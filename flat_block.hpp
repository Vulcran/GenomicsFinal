#ifndef FLAT_BLOCK_HPP
#define FLAT_BLOCK_HPP

#include <cstdint>
#include <cstddef>
#include <cstring>
#include <vector>

namespace flat_index {

/**
 * 64-byte block layout:
 * offset  size  field
 * 0       2    unitig_id            (uint16)
 * 2       1    dna_len_kmers        (uint8)
 * 3       1    n_utab               (uint8)
 * 4       1    n_etab               (uint8)
 * 5       1    flags                (uint8)
 * 6       2    next_block_offset/64 (uint16)
 * 8       D    DNA, 2-bit packed    (variable size)
 * 8+D    8*U   UTAB rows            (8 bytes/row)
 * ...    4*E   ETAB rows            (4 bytes/row)
 * ...   pad    zero padding to 64 B
 */

constexpr std::size_t BLOCK_SIZE = 64;
constexpr std::size_t HEADER_SIZE = 8;
constexpr std::size_t DNA_REGION_START = 8;
constexpr std::size_t DNA_REGION_SIZE = 40;

struct BlockHeader {
    uint16_t unitig_id;
    uint8_t  dna_len_kmers;
    uint8_t  n_utab;
    uint8_t  n_etab;
    uint8_t  flags;
    uint16_t next_block_offset_div_64;

    void pack(std::byte* dst) const {
        std::memcpy(dst + 0, &unitig_id, 2);
        std::memcpy(dst + 2, &dna_len_kmers, 1);
        std::memcpy(dst + 3, &n_utab, 1);
        std::memcpy(dst + 4, &n_etab, 1);
        std::memcpy(dst + 5, &flags, 1);
        std::memcpy(dst + 6, &next_block_offset_div_64, 2);
    }

    static BlockHeader unpack(const std::byte* src) {
        BlockHeader h;
        std::memcpy(&h.unitig_id, src + 0, 2);
        std::memcpy(&h.dna_len_kmers, src + 2, 1);
        std::memcpy(&h.n_utab, src + 3, 1);
        std::memcpy(&h.n_etab, src + 4, 1);
        std::memcpy(&h.flags, src + 5, 1);
        std::memcpy(&h.next_block_offset_div_64, src + 6, 2);
        return h;
    }
};

/**
 * UTAB row (8 bytes):
 * ref_id    : 16 bits
 * ref_pos   : 32 bits
 * entry_off : 8 bits
 * walk_len  : 7 bits
 * orient    : 1 bit (+ is 0, - is 1)
 */
struct UTabRow {
    uint16_t ref_id;
    uint32_t ref_pos;
    uint8_t  entry_off;
    uint8_t  walk_len; // 7 bits
    bool     orient;   // 1 bit

    void pack(std::byte* dst) const {
        std::memcpy(dst, &ref_id, 2);
        std::memcpy(dst + 2, &ref_pos, 4);
        dst[6] = static_cast<std::byte>(entry_off);
        dst[7] = static_cast<std::byte>((walk_len & 0x7F) | (orient ? 0x80 : 0x00));
    }

    static UTabRow unpack(const std::byte* src) {
        UTabRow r;
        std::memcpy(&r.ref_id, src, 2);
        std::memcpy(&r.ref_pos, src + 2, 4);
        r.entry_off = static_cast<uint8_t>(src[6]);
        r.walk_len = static_cast<uint8_t>(src[7]) & 0x7F;
        r.orient = (static_cast<uint8_t>(src[7]) & 0x80) != 0;
        return r;
    }
};

/**
 * ETAB row (4 bytes):
 * to_block_offset / 64 : 24 bits
 * ext_base             : 2 bits (A:0, C:1, G:2, T:3)
 * from_orient          : 1 bit (+:0, -:1)
 * to_orient            : 1 bit (+:0, -:1)
 * reserved             : 4 bits
 */
struct ETabRow {
    uint32_t to_block_offset_div_64; // 24 bits
    uint8_t  ext_base;               // 2 bits
    bool     from_orient;            // 1 bit
    bool     to_orient;              // 1 bit

    void pack(std::byte* dst) const {
        uint32_t packed = (to_block_offset_div_64 & 0xFFFFFF);
        packed |= (static_cast<uint32_t>(ext_base & 0x03) << 24);
        packed |= (from_orient ? 1u : 0u) << 26;
        packed |= (to_orient ? 1u : 0u) << 27;
        std::memcpy(dst, &packed, 4);
    }

    static ETabRow unpack(const std::byte* src) {
        ETabRow r;
        uint32_t packed;
        std::memcpy(&packed, src, 4);
        r.to_block_offset_div_64 = packed & 0xFFFFFF;
        r.ext_base = (packed >> 24) & 0x03;
        r.from_orient = (packed >> 26) & 0x01;
        r.to_orient = (packed >> 27) & 0x01;
        return r;
    }
};

} // namespace flat_index

#endif // FLAT_BLOCK_HPP
