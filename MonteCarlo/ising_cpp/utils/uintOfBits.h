#pragma once
#include <cstdint>

// Create a template that maps bit-widths to unsigned types
// clang-format off
template <unsigned Bytes> struct smallest_uint {
    static_assert(Bytes <= 16, "Bit width too large");
    using type = typename std::conditional_t<(Bytes <= 1), u_int8_t,
                 typename std::conditional_t<(Bytes <= 2), u_int16_t,
                 typename std::conditional_t<(Bytes <= 4), u_int32_t, 
                 typename std::conditional_t<(Bytes <= 8), u_int64_t, 
                 uint4>>>>;
};
// Alias for easier usage
template <std::size_t Bits> 
using uint_of_bits = typename smallest_uint<Bits>::type;
// clang-format on