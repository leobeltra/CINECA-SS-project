
#include "lut.h"

template <int NumNeighbours, typename T = double>
constexpr std::array<double, NumNeighbours> getProbabilityLUT(T J, T B_field, T beta) {
    std::array<double, NumNeighbours> lut;
    std::size_t lut_idx = 0;

    // All values of current_spin
    for (double current_spin = -1; current_spin <= 1; current_spin += 2) {

        // All values of neighbor_sum
        for (double sum = -NumNeighbours; sum < NumNeighbours; sum += 2) {
            // Calculate the change in energy (dE)
            T dE = T{2} * J * current_spin * sum + T{2} * B_field * current_spin;
            // Apply Metropolis algorithm
            lut[lut_idx] = exp(-beta * dE); // acceptance_prob
            lut_idx++;
        }
    }
    return lut;
}