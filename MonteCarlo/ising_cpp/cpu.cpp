#include "cpu.h"
#include <random>
#include <cmath>

void simulate_ising_cpu(real_t *spins, real_t *cpu_result, int N, int equil_steps, int M_sweep, 
                        real_t beta, real_t J, real_t B_field, real_t *cpu_acceptancy) {
    // Copy initial state to result
    for (int i = 0; i < N * N; i++) {
        cpu_result[i] = spins[i];
    }
    
    // Set up random number generation
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> pos_dist(0, N - 1);
    std::uniform_real_distribution<real_t> prob_dist(0.0, 1.0);
    
    int total_accepts = 0;
    int total_attempts = 0;
    
    // Equilibration phase
    for (int iter = 0; iter < equil_steps; iter++) {
        // Sweep the lattice
        for (int sweep = 0; sweep < N * N; sweep++) {
            // Select a random site
            int i = pos_dist(gen);
            int j = pos_dist(gen);
            int idx = i * N + j;
            
            // Get neighboring spins with periodic boundary conditions
            int up = ((i - 1 + N) % N) * N + j;
            int down = ((i + 1) % N) * N + j;
            int left = i * N + ((j - 1 + N) % N);
            int right = i * N + ((j + 1) % N);
            
            // Calculate energy change from flipping this spin
            real_t current_spin = cpu_result[idx];
            real_t neighbor_sum = cpu_result[up] + cpu_result[down] + cpu_result[left] + cpu_result[right];
            real_t energy_change = 2.0 * J * current_spin * neighbor_sum + 2.0 * B_field * current_spin;
            
            total_attempts++;
            
            // Metropolis acceptance criterion
            if (energy_change <= 0 || prob_dist(gen) < exp(-beta * energy_change)) {
                // Accept the flip
                cpu_result[idx] = -current_spin;
                total_accepts++;
            }
        }
    }
    
    // Measurement phase
    for (int iter = 0; iter < M_sweep; iter++) {
        // Sweep the lattice
        for (int sweep = 0; sweep < N * N; sweep++) {
            // Select a random site
            int i = pos_dist(gen);
            int j = pos_dist(gen);
            int idx = i * N + j;
            
            // Get neighboring spins with periodic boundary conditions
            int up = ((i - 1 + N) % N) * N + j;
            int down = ((i + 1) % N) * N + j;
            int left = i * N + ((j - 1 + N) % N);
            int right = i * N + ((j + 1) % N);
            
            // Calculate energy change from flipping this spin
            real_t current_spin = cpu_result[idx];
            real_t neighbor_sum = cpu_result[up] + cpu_result[down] + cpu_result[left] + cpu_result[right];
            real_t energy_change = 2.0 * J * current_spin * neighbor_sum + 2.0 * B_field * current_spin;
            
            total_attempts++;
            
            // Metropolis acceptance criterion
            if (energy_change <= 0 || prob_dist(gen) < exp(-beta * energy_change)) {
                // Accept the flip
                cpu_result[idx] = -current_spin;
                total_accepts++;
            }
        }
    }
    
    // Calculate and store acceptancy rate
    *cpu_acceptancy = static_cast<real_t>(total_accepts) / total_attempts;
}