#ifndef DNA_PACK_HPP
#define DNA_PACK_HPP

#include <cstdint>
#include <cstddef>
#include <string_view>
#include "kmer.hpp"

namespace flat_index {

/**
 * Packs bases into dst starting from the first bit.
 * Each base takes 2 bits.
 */
inline void pack_dna(std::string_view bases, std::byte* dst) {
    for (std::size_t i = 0; i < bases.size(); ++i) {
        uint8_t bits = base_to_bits(bases[i]);
        std::size_t byte_idx = (i * 2) / 8;
        std::size_t bit_idx = (i * 2) % 8;
        
        if (bit_idx == 0) dst[byte_idx] = std::byte{0};
        
        dst[byte_idx] |= static_cast<std::byte>(bits << (6 - bit_idx));
    }
}

/**
 * Extracts a 2-bit-packed k-mer starting at bit_off bits into src.
 * bit_off is relative to src.
 */
inline uint64_t read_kmer(const std::byte* src, std::size_t bit_off, int k) {
    uint64_t res = 0;
    for (int i = 0; i < k; ++i) {
        std::size_t current_bit = bit_off + i * 2;
        std::size_t byte_idx = current_bit / 8;
        std::size_t bit_in_byte = current_bit % 8;
        
        uint8_t byte_val = static_cast<uint8_t>(src[byte_idx]);
        uint8_t bits = (byte_val >> (6 - bit_in_byte)) & 0x03;
        res = (res << 2) | bits;
    }
    return res;
}

} // namespace flat_index

#endif // DNA_PACK_HPP
