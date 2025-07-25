#include "gpu.h"
#include <cuda_runtime.h>
#if DEBUG
#include <iostream>
#include <stdio.h>
#endif

// #define     CUDA_SAFE_CALL(call)                                                                       \
//     {                                                                                              \
//         call;                                                                                      \
//         cudaDeviceSynchronize();                                                                   \
//         cudaError_t err = cudaGetLastError();                                                      \
//         if (err != cudaSuccess) {                                                                  \
//             fprintf(stderr, "CUDA error in %s at %s:%d: %s\n", #call, __FILE__, __LINE__,          \
//                     cudaGetErrorString(err));                                                      \
//             exit(EXIT_FAILURE);                                                                    \
//         }                                                                                          \
//     }

#ifndef FAKE_CURAND
#include <curand_kernel.h>
#else
// Fake CURAND implementation for testing without dependencies
// Simple replacement for curandState
struct curandState {
    unsigned long long seed;
    unsigned long long sequence;
    unsigned long long offset;
    unsigned int state[16]; // Internal state array
};

// Fake curand_init implementation
__device__ void curand_init(unsigned long seed, unsigned long long sequence,
                            unsigned long long offset, curandState *state) {
    state->seed = seed;
    state->sequence = sequence;
    state->offset = offset;

    // Simple initialization of state array
    for (int i = 0; i < 16; i++) {
        state->state[i] = seed + sequence + i;
    }
}

// Fake curand_uniform implementation
__device__ float curand_uniform(curandState *state) {
    // Very simple LCG random number generator
    state->state[0] = 1664525 * state->state[0] + 1013904223;
    // Return a value between 0 and 1
    return (state->state[0] & 0x7FFFFFFF) / float(0x7FFFFFFF);
}

// Double precision version
__device__ double curand_uniform_double(curandState *state) {
    state->state[0] = 1664525 * state->state[0] + 1013904223;
    return (state->state[0] & 0x7FFFFFFF) / double(0x7FFFFFFF);
}
#endif

#define BLOCK_SIZE 16 // Define block size for CUDA kernels

// CUDA kernel for initializing random states, TODO: use real_t
__global__ void setup_rand_kernel(curandState *state, int N, unsigned long seed) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int idy = blockIdx.y * blockDim.y + threadIdx.y;

    if (idx >= N || idy >= N)
        return;

    int index = idy * N + idx;
    curand_init(seed, index, 0, &state[index]);
}

// CUDA kernel for Ising model simulation (checkerboard pattern, even sites)
template <typename T>
__global__ void ising_kernel_even(T *spins, T *new_spins, int N, T beta, T J, T B_field,
                                  curandState *rand_states, unsigned long long *accepted) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int idy = blockIdx.y * blockDim.y + threadIdx.y;

    // Skip odd sites in checkerboard pattern
    if ((idx + idy) % 2 != 0)
        return;
    // Skip out of bounds
    if (idx >= N || idy >= N)
        return;

    int index = idy * N + idx;
    curandState local_state = rand_states[index];

    // Calculate the indices of neighboring spins with periodic boundary conditions
    int left = idy * N + ((idx - 1 + N) % N);
    int right = idy * N + ((idx + 1) % N);
    int up = ((idy - 1 + N) % N) * N + idx;
    int down = ((idy + 1) % N) * N + idx;

    // Calculate energy change if we flip this spin
    const T current_spin = spins[index];
    T neighbor_sum = spins[left] + spins[right] + spins[up] + spins[down];
    T dE = T{2} * J * current_spin * neighbor_sum + T{2} * B_field * current_spin;

    // Apply Metropolis algorithm
    T acceptance_prob = exp(-beta * dE);
    const bool accept = curand_uniform(&local_state) < acceptance_prob;

    if (accept) {
        new_spins[index] = -current_spin;
        atomicAdd(accepted, 1ULL);
        spins[index] = -current_spin; // Update the spin state
    } else {
        new_spins[index] = current_spin;
    }

    // Update random state
    rand_states[index] = local_state;
}

// CUDA kernel for Ising model simulation (checkerboard pattern, odd sites)
template <typename T>
__global__ void ising_kernel_odd(T *spins, T *new_spins, int N, T beta, T J, T B_field,
                                 curandState *rand_states, unsigned long long *accepted) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int idy = blockIdx.y * blockDim.y + threadIdx.y;

    // Skip even sites in checkerboard pattern
    if ((idx + idy) % 2 == 0)
        return;

    // Skip out of bounds
    if (idx >= N || idy >= N)
        return;

    int index = idy * N + idx;
    curandState local_state = rand_states[index];

    // Calculate the indices of neighboring spins with periodic boundary conditions
    int left = idy * N + ((idx - 1 + N) % N);
    int right = idy * N + ((idx + 1) % N);
    int up = ((idy - 1 + N) % N) * N + idx;
    int down = ((idy + 1) % N) * N + idx;

    // Calculate energy change if we flip this spin
    const T current_spin = spins[index];
    T neighbor_sum = spins[left] + spins[right] + spins[up] + spins[down];
    T dE = T{2} * J * current_spin * neighbor_sum + T{2} * B_field * current_spin;

    // Apply Metropolis algorithm
    T acceptance_prob = exp(-beta * dE);
    const bool accept = curand_uniform(&local_state) < acceptance_prob;

    if (accept) {
        new_spins[index] = -current_spin;
        spins[index] = -current_spin; // Update the spin state
        atomicAdd(accepted, 1ULL);
    } else {
        new_spins[index] = current_spin;
    }

    // Update random state
    rand_states[index] = local_state;
}

// CUDA kernel for computing partial sums for Ising Hamiltonian
template <typename T>
__global__ void compute_hamiltonian_kernel(const T *spins, T *partial_sums, int N, T J, T B_field) {
    extern __shared__ T shared_data[];

    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int idy = blockIdx.y * blockDim.y + threadIdx.y;
    int tid = threadIdx.y * blockDim.x + threadIdx.x;

    T h_b = 0.0;
    T h_s_rows = 0.0;
    T h_s_cols = 0.0;

    if (idx < N && idy < N) {
        int index = idy * N + idx;

        // Calculate magnetic field contribution
        h_b = spins[index];

        // Calculate horizontal spin interactions
        h_s_rows = spins[index] * spins[idy * N + ((idx + 1) % N)];

        // Calculate vertical spin interactions
        h_s_cols = spins[index] * spins[((idy + 1) % N) * N + idx];
    }

    // Store partial results in shared memory
    shared_data[tid * 3 + 0] = h_b;
    shared_data[tid * 3 + 1] = h_s_rows;
    shared_data[tid * 3 + 2] = h_s_cols;
    __syncthreads();

    // Parallel reduction within block
    for (int s = blockDim.x * blockDim.y / 2; s > 0; s >>= 1) {
        if (tid < s) {
            shared_data[tid * 3 + 0] += shared_data[(tid + s) * 3 + 0];
            shared_data[tid * 3 + 1] += shared_data[(tid + s) * 3 + 1];
            shared_data[tid * 3 + 2] += shared_data[(tid + s) * 3 + 2];
        }
        __syncthreads();
    }

    // Write block results to global memory
    if (tid == 0) {
        int blockId = blockIdx.y * gridDim.x + blockIdx.x;
        partial_sums[blockId * 3 + 0] = shared_data[0]; // H_B contribution
        partial_sums[blockId * 3 + 1] = shared_data[1]; // H_S_rows contribution
        partial_sums[blockId * 3 + 2] = shared_data[2]; // H_S_cols contribution
    }
}

void simulate_ising_gpu(real_t *spins, real_t *result, int N, int equil_steps, int M_sweep,
                        real_t beta, real_t J, real_t B_field, real_t *acceptancy) {
    // Device memory allocation
    real_t *d_spins, *d_new_spins;
    size_t size = N * N * sizeof(real_t);
    cudaMalloc(&d_spins, size);
    cudaMalloc(&d_new_spins, size);

    // Copy input spins to device
    cudaMemcpy(d_spins, spins, size, cudaMemcpyHostToDevice);
    cudaMemcpy(d_new_spins, spins, size, cudaMemcpyHostToDevice);

    // Allocate memory for random states
    curandState *d_rand_states;
    cudaMalloc(&d_rand_states, N * N * sizeof(curandState));

    // Allocate memory for acceptancy counter
    unsigned long long *d_accepted, h_accepted = 0;
    cudaMalloc(&d_accepted, sizeof(unsigned long long));
    real_t temp_acceptancy = 0.0;

    // Setup grid and blocks for CUDA kernels
    dim3 block_size(BLOCK_SIZE, BLOCK_SIZE, 1);
    dim3 grid_size((N + block_size.x - 1) / block_size.x, (N + block_size.y - 1) / block_size.y, 1);

    int num_blocks = grid_size.x * grid_size.y;
    size_t shared_mem_size = block_size.x * block_size.y * 3 * sizeof(real_t);

    // Allocate memory for partial sums
    real_t *d_partial_sums;
    cudaMalloc(&d_partial_sums, num_blocks * 3 * sizeof(real_t));
    real_t *h_partial_sums;
    cudaMallocHost(&h_partial_sums, num_blocks * 3 * sizeof(real_t));

    // Initialize random generators
    setup_rand_kernel<<<grid_size, block_size>>>(d_rand_states, N, time(NULL));

#if DEBUG
    // print condition
    auto printCondition = [](int iter) { return (iter % 1 == 0 && iter < 5); };
#endif

    // Equilibration steps
    for (int iter = 0; iter < equil_steps; iter++) {
        cudaMemset(d_accepted, 0, sizeof(unsigned long long));
        ising_kernel_even<real_t><<<grid_size, block_size>>>(d_spins, d_new_spins, N, beta, J,
                                                             B_field, d_rand_states, d_accepted);
#if DEBUG
        cudaDeviceSynchronize();
        // Check cuda errors
        cudaError_t err = cudaGetLastError();
        if (err != cudaSuccess) {
            std::cerr << "CUDA error in ising_kernel_even: " << cudaGetErrorString(err)
                      << std::endl;
            exit(EXIT_FAILURE);
        }
#endif

        ising_kernel_odd<real_t><<<grid_size, block_size>>>(d_new_spins, d_spins, N, beta, J,
                                                            B_field, d_rand_states, d_accepted);

#if DEBUG
        // Copy d_spins to host to print the first 4x4 tile
        if (printCondition(iter)) { // Print every 10 iterations and first 500
            real_t *h_spins = new real_t[N * N];
            cudaMemcpy(h_spins, d_spins, size, cudaMemcpyDeviceToHost);
            std::cout << "Equilibration step " << iter << ": Spins (first 4x4 tile):" << std::endl;
            for (int i = 0; i < 4; i++) {
                std::cout << "\t";
                for (int j = 0; j < 4; j++) {
                    std::cout << h_spins[i * N + j] << " ";
                }
                std::cout << std::endl;
            }
            delete[] h_spins;
        }
#endif

        // Compute H_issing_2D
        // Launch kernel to compute partial sums
        compute_hamiltonian_kernel<real_t>
            <<<grid_size, block_size, shared_mem_size>>>(d_spins, d_partial_sums, N, J, B_field);

        // Copy partial sums back to host, automatic synchronization
        cudaMemcpy(h_partial_sums, d_partial_sums, num_blocks * 3 * sizeof(real_t),
                   cudaMemcpyDeviceToHost);

        // Sum up on host (could be optimized with another kernel for large grids)
        real_t h_b_sum = 0.0, h_s_rows_sum = 0.0, h_s_cols_sum = 0.0;
        for (int i = 0; i < num_blocks; i++) {
            h_b_sum += h_partial_sums[i * 3 + 0];
            h_s_rows_sum += h_partial_sums[i * 3 + 1];
            h_s_cols_sum += h_partial_sums[i * 3 + 2];
        }

        // Calculate final Hamiltonian
        real_t h_ising = -B_field * h_b_sum - J * (h_s_rows_sum + h_s_cols_sum);

#if DEBUG
        if (printCondition(iter)) {
            // Print Hamiltonian
            std::cout << "Equilibration step " << iter << ": H_ising = " << h_ising << std::endl;
            // Print acceptance count
            cudaMemcpy(&h_accepted, d_accepted, sizeof(unsigned long long), cudaMemcpyDeviceToHost);
            cudaDeviceSynchronize();
            std::cout << "Acceptance count: " << h_accepted << " N*N/2=" << N * N / 2 << std::endl;
        }
#endif
    }

    // Copy acceptance count to host
    cudaMemcpy(&h_accepted, d_accepted, sizeof(unsigned long long), cudaMemcpyDeviceToHost);
    temp_acceptancy = (real_t)h_accepted / (real_t)(N * N);

#if DEBUG
    // Print acceptance rate after equilibration
    std::cout << "Acceptance rate after equilibration: " << temp_acceptancy << std::endl;
#endif

    // Measurement steps
    for (int iter = 0; iter < M_sweep; iter++) {

        // Reset acceptance counter for measurement phase
        cudaMemset(d_accepted, 0, sizeof(unsigned long long));

        ising_kernel_even<real_t><<<grid_size, block_size>>>(d_spins, d_new_spins, N, beta, J,
                                                             B_field, d_rand_states, d_accepted);
        ising_kernel_odd<real_t><<<grid_size, block_size>>>(d_new_spins, d_spins, N, beta, J,
                                                            B_field, d_rand_states, d_accepted);

        // Copy acceptance count to host, implicit synchronization
        cudaMemcpy(&h_accepted, d_accepted, sizeof(unsigned long long), cudaMemcpyDeviceToHost);
        temp_acceptancy = (real_t)h_accepted / (real_t)(N * N);
        *acceptancy = *acceptancy + temp_acceptancy;
    }

    // Copy results back to host
    cudaMemcpy(result, d_spins, size, cudaMemcpyDeviceToHost);
    cudaMemcpy(&h_accepted, d_accepted, sizeof(unsigned long long), cudaMemcpyDeviceToHost);

    // Calculate acceptancy rate
    *acceptancy = *acceptancy / (real_t)(M_sweep);

    // Free device memory
    cudaFree(d_spins);
    cudaFree(d_new_spins);
    cudaFree(d_rand_states);
    cudaFree(d_accepted);
    cudaFree(d_partial_sums);
    // Free host memory
    cudaFreeHost(h_partial_sums);
}