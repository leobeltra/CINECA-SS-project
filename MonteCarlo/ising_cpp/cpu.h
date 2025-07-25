#pragma once
#include "general.h"

void simulate_ising_cpu(real_t *spins, real_t *cpu_result, int N, int equil_steps, int M_sweep, 
                        real_t beta, real_t J, real_t B_field, real_t *cpu_acceptancy);