
#pragma once
#include <array>

template <int NumNeighbours, typename T = double>
constexpr std::array<double, 2 * NumNeighbours> getProbabilityLUT(T J, T B_field, T beta);