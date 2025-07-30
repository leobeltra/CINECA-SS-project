#include "../gpu.h"
#include "../utils/lut.h"
#include "../utils/uintOfBits.h"
#include <cuda_runtime.h>
#if DEBUG
#include <iostream>
#include <stdio.h>
#endif
#include <curand_kernel.h>

#define BLOCK_SIZE 16 // Define block size for CUDA kernels

#define NUM_PROBABILITIES 8 // 2 (s_i combinations) * 4 (sum s_j combinations)

#define FLIP_WITH_PROBABILITY(el_yx)                                                               \
    probability = probabilities[spinToUpdate * NUM_PROBABILITIES / 2 + neighborsSum];              \
    accept = curand_uniform(&generatorState) < probability;                                        \
    if (accept) {                                                                                  \
        el_yx = !spinToUpdate; /* Flip the spin */                                                 \
        localAccepted++;                                                                           \
    }

struct TileEvenElements {
    union {
        u_int64_t tile;
        struct {
            u_int8_t el_00; // el_yx
            u_int8_t el_02;
            u_int8_t el_11;
            u_int8_t el_13;
            u_int8_t el_20;
            u_int8_t el_22;
            u_int8_t el_31;
            u_int8_t el_33;
        };
    };
};

struct TileOddElements {
    union {
        u_int64_t tile;
        struct {
            u_int8_t el_01; // el_yx
            u_int8_t el_03;
            u_int8_t el_10;
            u_int8_t el_12;
            u_int8_t el_21;
            u_int8_t el_23;
            u_int8_t el_30;
            u_int8_t el_32;
        };
    };
};

// CUDA kernel for initializing random states, TODO: use real_t
__global__ void setup_rand_kernel(curandState *state, int N, unsigned long seed) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int idy = blockIdx.y * blockDim.y + threadIdx.y;

    if (idx >= N || idy >= N)
        return;

    int index = idy * N + idx;
    curand_init(seed, index, 0, &state[index]);
}

struct prob {
    real_t data[NUM_PROBABILITIES];
};

// CUDA kernel for Ising model simulation (checkerboard pattern, update even sites)
template <typename T>
__global__ void ising_kernel_even(TileEvenElements *even_spins,
                                  const TileOddElements *odd_spins,
                                  int qN,                                  // grid size divided by 4
                                  const __grid_constant__ prob probStruct, // pass lut by value
                                  curandState *rand_states,
                                  unsigned long long *accepted) {
    const int idxX = blockIdx.x * blockDim.x + threadIdx.x;
    const int idxY = blockIdx.y * blockDim.y + threadIdx.y;

    // Skip out of bounds
    if (idxX >= qN || idxY >= qN)
        return;

    const int linIdx = idxY * qN + idxX;
    curandState generatorState = rand_states[linIdx];
    unsigned int localAccepted = 0;

    // Calculate the indices of neighboring tiles with periodic boundary conditions
    const int leftIdx = idxY * qN + ((idxX - 1 + qN) % qN);
    const int rightIdx = idxY * qN + ((idxX + 1) % qN);
    const int topIdx = ((idxY - 1 + qN) % qN) * qN + idxX;
    const int downIdx = ((idxY + 1) % qN) * qN + idxX;
    const auto &probabilities = probStruct.data; // access the LUT data

    // ---------------------------------------------
    // Update spins of all elements in the even tile
    // ---------------------------------------------

    // Load old spinds to update
    TileEvenElements localTileEven = even_spins[linIdx]; // elements are 0 for -1, 1 for +1

    // start updating the two central elements: el_11 and el_22
    // load localTileOdd elements
    TileOddElements localTileOdd = odd_spins[linIdx];
    // el_11
    u_int16_t spinToUpdate = localTileEven.el_11;
    u_int16_t neighborsSum =
        localTileOdd.el_10 + localTileOdd.el_12 + localTileOdd.el_21 + localTileOdd.el_01;
    // this repeating pattern will be encapsulated in a macro: FLIP_WITH_PROBABILITY(el_xy)
    T probability = probabilities[spinToUpdate * NUM_PROBABILITIES / 2 + neighborsSum];
    bool accept = curand_uniform(&generatorState) < probability;
    if (accept) {
        localTileEven.el_11 = !spinToUpdate; // flip the spin
        localAccepted++;
    }
    // el_22
    spinToUpdate = localTileEven.el_22;
    neighborsSum =
        localTileOdd.el_21 + localTileOdd.el_23 + localTileOdd.el_32 + localTileOdd.el_12;
    FLIP_WITH_PROBABILITY(localTileEven.el_22);

    // update the elements on the left: el_00, el_20
    // el_20
    // need to load the left odd tile
    TileOddElements leftTileOdd = odd_spins[leftIdx];
    spinToUpdate = localTileEven.el_20;
    neighborsSum = leftTileOdd.el_23 + localTileOdd.el_21 + localTileOdd.el_30 + localTileOdd.el_10;
    FLIP_WITH_PROBABILITY(localTileEven.el_20);
    // el_00
    // need to have the left and top odd tiles
    TileOddElements upTileOdd = odd_spins[topIdx];
    spinToUpdate = localTileEven.el_00;
    neighborsSum = leftTileOdd.el_03 + localTileOdd.el_01 + localTileOdd.el_10 + upTileOdd.el_30;
    FLIP_WITH_PROBABILITY(localTileEven.el_00);

    // update the elements on the top: el_02
    // el_02
    spinToUpdate = localTileEven.el_02;
    neighborsSum = localTileOdd.el_01 + localTileOdd.el_03 + upTileOdd.el_32 + localTileOdd.el_12;
    FLIP_WITH_PROBABILITY(localTileEven.el_02);

    // update the elements on the right: el_13, el_33
    // el_13
    // need to load the right odd tile
    TileOddElements rightTileOdd = odd_spins[rightIdx];
    spinToUpdate = localTileEven.el_13;
    neighborsSum =
        localTileOdd.el_12 + rightTileOdd.el_10 + localTileOdd.el_03 + localTileOdd.el_32;
    FLIP_WITH_PROBABILITY(localTileEven.el_13);
    // el_33
    // need to have the right and down odd tiles
    TileOddElements downTileOdd = odd_spins[downIdx];
    spinToUpdate = localTileEven.el_33;
    neighborsSum = localTileOdd.el_32 + rightTileOdd.el_30 + localTileOdd.el_23 + downTileOdd.el_03;
    FLIP_WITH_PROBABILITY(localTileEven.el_33);

    // update the element on the bottom: el_31
    // el_31
    spinToUpdate = localTileEven.el_31;
    neighborsSum = localTileOdd.el_30 + localTileOdd.el_32 + localTileOdd.el_21 + downTileOdd.el_01;
    FLIP_WITH_PROBABILITY(localTileEven.el_31);

    // --------------------------------
    // Update global variables memory
    // --------------------------------

    // Write back the updated tile
    even_spins[linIdx] = localTileEven;

    // Update the acceptance count
    atomicAdd(accepted, localAccepted);

    // Update random state
    rand_states[linIdx] = generatorState;
}

// CUDA kernel for Ising model simulation (checkerboard pattern, update odd sites)
// [generated by copilot based on the other]
template <typename T>
__global__ void ising_kernel_odd(TileOddElements *odd_spins,
                                 const TileEvenElements *even_spins,
                                 int qN,                                  // grid size divided by 4
                                 const __grid_constant__ prob probStruct, // passed by value
                                 curandState *rand_states,
                                 unsigned long long *accepted) {
    const int idxX = blockIdx.x * blockDim.x + threadIdx.x;
    const int idxY = blockIdx.y * blockDim.y + threadIdx.y;

    // Skip out of bounds
    if (idxX >= qN || idxY >= qN)
        return;

    const int linIdx = idxY * qN + idxX;
    curandState generatorState = rand_states[linIdx];
    unsigned int localAccepted = 0;
    T probability;
    bool accept;
    const auto &probabilities = probStruct.data; // access the LUT data

    // Calculate the indices of neighboring tiles with periodic boundary conditions
    const int leftIdx = idxY * qN + ((idxX - 1 + qN) % qN);
    const int rightIdx = idxY * qN + ((idxX + 1) % qN);
    const int topIdx = ((idxY - 1 + qN) % qN) * qN + idxX;
    const int downIdx = ((idxY + 1) % qN) * qN + idxX;

    // ---------------------------------------------
    // Update spins of all elements in the odd tile
    // ---------------------------------------------

    // Load old spins to update
    TileOddElements localTileOdd = odd_spins[linIdx];

    // Load local and neighboring even elements
    TileEvenElements localTileEven = even_spins[linIdx];
    TileEvenElements leftTileEven = even_spins[leftIdx];
    TileEvenElements topTileEven = even_spins[topIdx];
    TileEvenElements rightTileEven = even_spins[rightIdx];
    TileEvenElements downTileEven = even_spins[downIdx];

    // Start updating the central elements: el_10, el_12, el_21, el_23
    // el_10
    u_int16_t spinToUpdate = localTileOdd.el_10;
    u_int16_t neighborsSum =
        localTileEven.el_00 + localTileEven.el_20 + localTileEven.el_11 + leftTileEven.el_13;
    FLIP_WITH_PROBABILITY(localTileOdd.el_10);

    // el_12
    spinToUpdate = localTileOdd.el_12;
    neighborsSum =
        localTileEven.el_02 + localTileEven.el_22 + localTileEven.el_11 + localTileEven.el_13;
    FLIP_WITH_PROBABILITY(localTileOdd.el_12);

    // el_21
    spinToUpdate = localTileOdd.el_21;
    neighborsSum =
        localTileEven.el_11 + localTileEven.el_31 + localTileEven.el_20 + localTileEven.el_22;
    FLIP_WITH_PROBABILITY(localTileOdd.el_21);

    // el_23
    spinToUpdate = localTileOdd.el_23;
    neighborsSum =
        localTileEven.el_13 + localTileEven.el_33 + localTileEven.el_22 + rightTileEven.el_20;
    FLIP_WITH_PROBABILITY(localTileOdd.el_23);

    // Update the edge elements: el_01, el_03, el_30, el_32
    // el_01
    spinToUpdate = localTileOdd.el_01;
    neighborsSum =
        localTileEven.el_00 + localTileEven.el_02 + localTileEven.el_11 + topTileEven.el_31;
    FLIP_WITH_PROBABILITY(localTileOdd.el_01);

    // el_03
    spinToUpdate = localTileOdd.el_03;
    neighborsSum =
        localTileEven.el_02 + rightTileEven.el_00 + localTileEven.el_13 + topTileEven.el_33;
    FLIP_WITH_PROBABILITY(localTileOdd.el_03);

    // el_30
    spinToUpdate = localTileOdd.el_30;
    neighborsSum =
        localTileEven.el_20 + localTileEven.el_31 + downTileEven.el_00 + leftTileEven.el_33;
    FLIP_WITH_PROBABILITY(localTileOdd.el_30);

    // el_32
    spinToUpdate = localTileOdd.el_32;
    neighborsSum =
        localTileEven.el_22 + localTileEven.el_33 + downTileEven.el_02 + localTileEven.el_31;
    FLIP_WITH_PROBABILITY(localTileOdd.el_32);

    // --------------------------------
    // Update global variables memory
    // --------------------------------

    // Write back the updated tile
    odd_spins[linIdx] = localTileOdd;

    // Update the acceptance count
    atomicAdd(accepted, localAccepted);

    // Update random state
    rand_states[linIdx] = generatorState;
}
#undef probababilities

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

// CUDA kernel to split the flat array of spins into even and odd tiles
// [generated by copilot sorry for this mess, but should work]
__global__ void copySpinsToDevice(const real_t *spins,
                                  int N,
                                  TileEvenElements *even_spins,
                                  TileOddElements *odd_spins) {
    // Calculate the tile index this thread is responsible for
    int tileX = blockIdx.x * blockDim.x + threadIdx.x;
    int tileY = blockIdx.y * blockDim.y + threadIdx.y;

    // Check if this thread is within bounds
    if (tileX >= N / 4 || tileY >= N / 4)
        return;

    // Calculate the linear index of this tile
    int tileIdx = tileY * (N / 4) + tileX;

    // Calculate the top-left corner of this tile in the original grid
    int startX = tileX * 4;
    int startY = tileY * 4;

    // Initialize the even and odd elements structures
    TileEvenElements evnTile;
    TileOddElements oddTile;

    // Fill the elements
    evnTile.el_00 = spins[(startY + 0) * N + (startX + 0)] > 0 ? 1 : 0;
    oddTile.el_01 = spins[(startY + 0) * N + (startX + 1)] > 0 ? 1 : 0;
    evnTile.el_02 = spins[(startY + 0) * N + (startX + 2)] > 0 ? 1 : 0;
    oddTile.el_03 = spins[(startY + 0) * N + (startX + 3)] > 0 ? 1 : 0;

    oddTile.el_10 = spins[(startY + 1) * N + (startX + 0)] > 0 ? 1 : 0;
    evnTile.el_11 = spins[(startY + 1) * N + (startX + 1)] > 0 ? 1 : 0;
    oddTile.el_12 = spins[(startY + 1) * N + (startX + 2)] > 0 ? 1 : 0;
    evnTile.el_13 = spins[(startY + 1) * N + (startX + 3)] > 0 ? 1 : 0;

    evnTile.el_20 = spins[(startY + 2) * N + (startX + 0)] > 0 ? 1 : 0;
    oddTile.el_21 = spins[(startY + 2) * N + (startX + 1)] > 0 ? 1 : 0;
    evnTile.el_22 = spins[(startY + 2) * N + (startX + 2)] > 0 ? 1 : 0;
    oddTile.el_23 = spins[(startY + 2) * N + (startX + 3)] > 0 ? 1 : 0;

    oddTile.el_30 = spins[(startY + 3) * N + (startX + 0)] > 0 ? 1 : 0;
    evnTile.el_31 = spins[(startY + 3) * N + (startX + 1)] > 0 ? 1 : 0;
    evnTile.el_33 = spins[(startY + 3) * N + (startX + 3)] > 0 ? 1 : 0;
    oddTile.el_32 = spins[(startY + 3) * N + (startX + 2)] > 0 ? 1 : 0;

    // Write the tiles back to global memory
    evnTile.tile = evnTile.tile; // added by me, check ptx if this is needed
    oddTile.tile = oddTile.tile;
    even_spins[tileIdx] = evnTile;
    odd_spins[tileIdx] = oddTile;
}

// CUDA kernel to place the even and odd tiles back to the original flat array of spins
// [generated by copilot]
__global__ void copySpinsToHost(const TileEvenElements *even_spins,
                                const TileOddElements *odd_spins,
                                int N,
                                real_t *spins) {
    // Calculate the tile index this thread is responsible for
    int tileX = blockIdx.x * blockDim.x + threadIdx.x;
    int tileY = blockIdx.y * blockDim.y + threadIdx.y;

    // Check if this thread is within bounds
    if (tileX >= N / 4 || tileY >= N / 4)
        return;

    // Calculate the linear index of this tile
    int tileIdx = tileY * (N / 4) + tileX;

    // Calculate the top-left corner of this tile in the original grid
    int startX = tileX * 4;
    int startY = tileY * 4;

    // Get the even and odd tile structures
    TileEvenElements evnTile = even_spins[tileIdx];
    TileOddElements oddTile = odd_spins[tileIdx];

    // Map the binary values back to -1.0/+1.0 real_t values and write to the spins array
    spins[(startY + 0) * N + (startX + 0)] = evnTile.el_00 ? 1.0 : -1.0;
    spins[(startY + 0) * N + (startX + 1)] = oddTile.el_01 ? 1.0 : -1.0;
    spins[(startY + 0) * N + (startX + 2)] = evnTile.el_02 ? 1.0 : -1.0;
    spins[(startY + 0) * N + (startX + 3)] = oddTile.el_03 ? 1.0 : -1.0;

    spins[(startY + 1) * N + (startX + 0)] = oddTile.el_10 ? 1.0 : -1.0;
    spins[(startY + 1) * N + (startX + 1)] = evnTile.el_11 ? 1.0 : -1.0;
    spins[(startY + 1) * N + (startX + 2)] = oddTile.el_12 ? 1.0 : -1.0;
    spins[(startY + 1) * N + (startX + 3)] = evnTile.el_13 ? 1.0 : -1.0;

    spins[(startY + 2) * N + (startX + 0)] = evnTile.el_20 ? 1.0 : -1.0;
    spins[(startY + 2) * N + (startX + 1)] = oddTile.el_21 ? 1.0 : -1.0;
    spins[(startY + 2) * N + (startX + 2)] = evnTile.el_22 ? 1.0 : -1.0;
    spins[(startY + 2) * N + (startX + 3)] = oddTile.el_23 ? 1.0 : -1.0;

    spins[(startY + 3) * N + (startX + 0)] = oddTile.el_30 ? 1.0 : -1.0;
    spins[(startY + 3) * N + (startX + 1)] = evnTile.el_31 ? 1.0 : -1.0;
    spins[(startY + 3) * N + (startX + 2)] = oddTile.el_32 ? 1.0 : -1.0;
    spins[(startY + 3) * N + (startX + 3)] = evnTile.el_33 ? 1.0 : -1.0;
}

void simulate_ising_gpu(real_t *spins,
                        real_t *result,
                        int N,
                        int equil_steps,
                        int M_sweep,
                        real_t beta,
                        real_t J,
                        real_t B_field,
                        real_t *acceptancy) {
    // Setup grid and blocks for CUDA kernels
    dim3 block_size(BLOCK_SIZE, BLOCK_SIZE, 1);
    dim3 grid_size((N / 4 + block_size.x - 1) / block_size.x,
                   (N / 4 + block_size.y - 1) / block_size.y, 1);

    // Device memory allocation
    TileEvenElements *even_spins;
    TileOddElements *odd_spins;
    static_assert(sizeof(TileEvenElements) == sizeof(TileOddElements), "Tile sizes must match");
    static_assert(sizeof(TileEvenElements) == 8,
                  "Don't forget to also change this size, a lot of places assume 4x4 tiles");
    size_t size = N / 4 * N / 4 * sizeof(TileEvenElements);
    cudaMalloc(&even_spins, size);
    cudaMalloc(&odd_spins, size);

    // Use a kernel to split spins to two tiles.
    // Don't thing it's needed, we can rand initial state inside the kernel
    real_t *d_spins;
    cudaMalloc(&d_spins, N * N * sizeof(real_t));
    cudaMemcpy(d_spins, spins, N * N * sizeof(real_t), cudaMemcpyHostToDevice);

    // could run on a different stream and synch after setup rand
    copySpinsToDevice<<<block_size, grid_size>>>(d_spins, N, even_spins, odd_spins);

    // Allocate memory for random states
    curandState *d_rand_states;
    cudaMalloc(&d_rand_states, N / 4 * N / 4 * sizeof(curandState));

    // Allocate memory for acceptancy counter
    unsigned long long *d_accepted, h_accepted = 0;
    cudaMalloc(&d_accepted, sizeof(unsigned long long));
    real_t temp_acceptancy = 0.0;

    int num_blocks = grid_size.x * grid_size.y;
    size_t shared_mem_size = block_size.x * block_size.y * 3 * sizeof(real_t);

    // Allocate memory for partial sums
    real_t *d_partial_sums;
    cudaMalloc(&d_partial_sums, num_blocks * 3 * sizeof(real_t));
    real_t *h_partial_sums;
    cudaMallocHost(&h_partial_sums, num_blocks * 3 * sizeof(real_t));

    // Initialize random generators
    setup_rand_kernel<<<grid_size, block_size>>>(d_rand_states, N / 4, time(NULL));

    int index = 0;
    prob probabilities;
    for (auto prob : getProbabilityLUT<NUM_PROBABILITIES, real_t>(J, B_field, beta)) {
        probabilities.data[index++] = prob;
    }

#if DEBUG
    // print condition
    auto printCondition = [](int iter) { return (iter % 1 == 0 && iter < 5); };
#endif
    // Equilibration steps
    for (int iter = 0; iter < equil_steps; iter++) {
        cudaMemset(d_accepted, 0, sizeof(unsigned long long));
        ising_kernel_even<real_t><<<grid_size, block_size>>>(
            even_spins, odd_spins, N / 4, probabilities, d_rand_states, d_accepted);
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

        ising_kernel_odd<real_t><<<grid_size, block_size>>>(
            odd_spins, even_spins, N / 4, probabilities, d_rand_states, d_accepted);

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

        // TODO: change this
        // // Compute H_issing_2D
        // // Launch kernel to compute partial sums
        // compute_hamiltonian_kernel<real_t>
        //     <<<grid_size, block_size, shared_mem_size>>>(d_spins, d_partial_sums, N, J, B_field);

        // // Copy partial sums back to host, automatic synchronization
        // cudaMemcpy(h_partial_sums, d_partial_sums, num_blocks * 3 * sizeof(real_t),
        //            cudaMemcpyDeviceToHost);

        // // Sum up on host (could be optimized with another kernel for large grids)
        // real_t h_b_sum = 0.0, h_s_rows_sum = 0.0, h_s_cols_sum = 0.0;
        // for (int i = 0; i < num_blocks; i++) {
        //     h_b_sum += h_partial_sums[i * 3 + 0];
        //     h_s_rows_sum += h_partial_sums[i * 3 + 1];
        //     h_s_cols_sum += h_partial_sums[i * 3 + 2];
        // }

        // // Calculate final Hamiltonian
        // real_t h_ising = -B_field * h_b_sum - J * (h_s_rows_sum + h_s_cols_sum);

#if DEBUG
        // if (printCondition(iter)) {
        //     // Print Hamiltonian
        //     std::cout << "Equilibration step " << iter << ": H_ising = " << h_ising << std::endl;
        //     // Print acceptance count
        //     cudaMemcpy(&h_accepted, d_accepted, sizeof(unsigned long long),
        //     cudaMemcpyDeviceToHost); cudaDeviceSynchronize(); std::cout << "Acceptance count: "
        //     << h_accepted << " N*N/2=" << N * N / 2 << std::endl;
        // }
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

        ising_kernel_even<real_t><<<grid_size, block_size>>>(
            even_spins, odd_spins, N / 4, probabilities, d_rand_states, d_accepted);
        ising_kernel_odd<real_t><<<grid_size, block_size>>>(
            odd_spins, even_spins, N / 4, probabilities, d_rand_states, d_accepted);

        // Copy acceptance count to host, implicit synchronization
        cudaMemcpy(&h_accepted, d_accepted, sizeof(unsigned long long), cudaMemcpyDeviceToHost);
        temp_acceptancy = (real_t)h_accepted / (real_t)(N * N);
        *acceptancy = *acceptancy + temp_acceptancy;
    }

    // Copy results back to host
    copySpinsToHost<<<grid_size, block_size>>>(even_spins, odd_spins, N, d_spins);
    cudaMemcpy(result, d_spins, size, cudaMemcpyDeviceToHost);
    cudaMemcpy(&h_accepted, d_accepted, sizeof(unsigned long long), cudaMemcpyDeviceToHost);

    // Calculate acceptancy rate
    *acceptancy = *acceptancy / (real_t)(M_sweep);

    // Free device memory
    cudaFree(even_spins);
    cudaFree(odd_spins);
    cudaFree(d_rand_states);
    cudaFree(d_accepted);
    cudaFree(d_partial_sums);
    cudaFree(d_spins);
    // Free host memory
    cudaFreeHost(h_partial_sums);
}