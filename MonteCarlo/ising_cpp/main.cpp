#include "cpu.h"
#include "general.h"
#include "gpu.h"
#include <chrono>
#include <cmath>
#include <cuda_runtime.h>
#include <iomanip>
#include <iostream>
#include <vector>

int main(int argc, char *argv[]) {
    // Physical parameters
    constexpr int N = 256;            // Lattice size (NxN)
    constexpr int equil_steps = 1000; // Equilibration steps
    constexpr int M_sweep = 10000;    // Measurement steps after equilibration
    constexpr real_t J = 0.1;         // Coupling constant
    constexpr real_t B_field = 0.0;   // External magnetic field
    constexpr real_t beta = 4.5;      // Inverse temperature (1/kT), default around critical point

    std::cout << "2D Ising Model Simulation" << std::endl;
#if USE_DOUBLE_PRECISION
    std::cout << "Using double precision" << std::endl;
#else
    std::cout << "Using single precision (float)" << std::endl;
#endif
    std::cout << "Lattice size: " << N << "x" << N << std::endl;
    std::cout << "Equilibration steps: " << equil_steps << std::endl;
    std::cout << "Measurement steps: " << M_sweep << std::endl;
    std::cout << "Beta (1/kT): " << beta << std::endl;
    std::cout << "Coupling constant (J): " << J << std::endl;
    std::cout << "External field (B): " << B_field << std::endl;

    // Use pinned memory for better CPU-GPU transfer performance
    real_t *spins = nullptr;
    real_t *cpu_result = nullptr;
    real_t *gpu_result = nullptr;
    real_t *cpu_acceptancy = nullptr;
    real_t *gpu_acceptancy = nullptr;

    // Allocate pinned memory
    cudaMallocHost(&spins, N * N * sizeof(real_t));
    cudaMallocHost(&cpu_result, N * N * sizeof(real_t));
    cudaMallocHost(&gpu_result, N * N * sizeof(real_t));
    cudaMallocHost(&cpu_acceptancy, sizeof(real_t));
    cudaMallocHost(&gpu_acceptancy, sizeof(real_t));

    double cpu_time_total = 0.0;
    double gpu_time_total = 0.0;

    // Run simulations multiple times to get average performance
    for (int run = 0; run < 1 + NUM_EXE_TO_MEASURE; run++) {
        std::cout << "Run " << run << (run == 0 ? " (warmup)" : "") << std::endl;

        // Initialize spins randomly
        for (int i = 0; i < N * N; i++) {
            spins[i] = (rand() % 2 == 0) ? 1.0 : -1.0;
        }

        // Time CPU implementation
#if RUN_CPU
        {
            auto start = std::chrono::high_resolution_clock::now();

            // Call CPU implementation with new parameters
            simulate_ising_cpu(spins, cpu_result, N, equil_steps, M_sweep, beta, J, B_field,
                               cpu_acceptancy);

            auto end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = end - start;

            if (run > 0) { // Skip first run when calculating average
                cpu_time_total += elapsed.count();
            }

            if (run == NUM_EXE_TO_MEASURE) { // Only print for the last run
                std::cout << std::fixed << std::setprecision(4);
                std::cout << "CPU acceptancy rate: " << *cpu_acceptancy << std::endl;
            }
        }
#endif
        // Time GPU implementation
        {
            auto start = std::chrono::high_resolution_clock::now();

            // Call GPU implementation with new parameters
            simulate_ising_gpu(spins, gpu_result, N, equil_steps, M_sweep, beta, J, B_field,
                               gpu_acceptancy);

            auto end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = end - start;

            if (run > 0) { // Skip first run when calculating average
                gpu_time_total += elapsed.count();
            }

            if (run == NUM_EXE_TO_MEASURE) { // Only print for the last run
                std::cout << "GPU acceptancy rate: " << *gpu_acceptancy << std::endl;
            }
        }
    }

    // Calculate average times
    double cpu_time = cpu_time_total / NUM_EXE_TO_MEASURE;
    double gpu_time = gpu_time_total / NUM_EXE_TO_MEASURE;

    std::cout << "\nPerformance results (averaged over " << NUM_EXE_TO_MEASURE
              << " runs):" << std::endl;
#if RUN_CPU
    std::cout << "CPU implementation took " << cpu_time << " seconds (avg)" << std::endl;
#endif
    std::cout << "GPU implementation took " << gpu_time << " seconds (avg)" << std::endl;

    // Calculate magnetization from final state
    double cpu_magnetization = 0.0;
    double gpu_magnetization = 0.0;

    for (int i = 0; i < N * N; i++) {
        cpu_magnetization += cpu_result[i];
        gpu_magnetization += gpu_result[i];
    }

    cpu_magnetization /= (N * N);
    gpu_magnetization /= (N * N);

#if RUN_CPU
    std::cout << "CPU final magnetization: " << cpu_magnetization << std::endl;
#else
    std::cout << "CPU final magnetization: N/A (CPU run disabled)" << std::endl;
#endif
    std::cout << "GPU final magnetization: " << gpu_magnetization << std::endl;

    // Calculate theoretical magnetization for comparison (if B_field is non-zero)
    if (B_field != 0.0) {
        double a = sinh(beta * B_field) +
                   (sinh(beta * B_field) * cosh(beta * B_field)) /
                       (sqrt(sinh(beta * B_field) * sinh(beta * B_field) + exp(-4 * beta * J)));
        double b = cosh(beta * B_field) +
                   sqrt(exp(-4 * beta * J) + sinh(beta * B_field) * sinh(beta * B_field));
        double magnetization_th = a / b;
        std::cout << "Theoretical magnetization: " << magnetization_th << std::endl;
    } else {
        std::cout << "Theoretical magnetization: 0.0 (no external field)" << std::endl;
    }

#if RUN_CPU
    std::cout << "Speedup: " << (cpu_time / gpu_time) << "x" << std::endl;
#else
    std::cout << "Speedup: N/A (CPU run disabled)" << std::endl;
#endif

    // Free pinned memory
    cudaFreeHost(spins);
    cudaFreeHost(cpu_result);
    cudaFreeHost(gpu_result);
    cudaFreeHost(cpu_acceptancy);
    cudaFreeHost(gpu_acceptancy);

    return 0;
}