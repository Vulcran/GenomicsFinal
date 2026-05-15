#ifndef KMER_HPP
#define KMER_HPP

#include <cstdint>
#include <string>
#include <algorithm>
#include <stdexcept>

namespace flat_index {

inline uint8_t base_to_bits(char base) {
    switch (base) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default: throw std::runtime_error("Invalid DNA base");
    }
}

inline char bits_to_base(uint8_t bits) {
    switch (bits & 0x03) {
        case 0: return 'A';
        case 1: return 'C';
        case 2: return 'G';
        case 3: return 'T';
        default: return 'N';
    }
}

inline uint64_t encode(const std::string& kmer) {
    uint64_t res = 0;
    for (char c : kmer) {
        res = (res << 2) | base_to_bits(c);
    }
    return res;
}

inline std::string decode(uint64_t kmer, int k) {
    std::string res(k, ' ');
    for (int i = 0; i < k; ++i) {
        res[k - 1 - i] = bits_to_base(kmer & 0x03);
        kmer >>= 2;
    }
    return res;
}

inline uint64_t reverse_complement(uint64_t kmer, int k) {
    uint64_t res = 0;
    for (int i = 0; i < k; ++i) {
        uint8_t bits = kmer & 0x03;
        uint8_t rc_bits = bits ^ 0x03; // 00->11 (A->T), 01->10 (C->G), etc.
        res = (res << 2) | rc_bits;
        kmer >>= 2;
    }
    return res;
}

inline uint64_t canonical(uint64_t kmer, int k) {
    return std::min(kmer, reverse_complement(kmer, k));
}

inline bool is_canonical(uint64_t kmer, int k) {
    return kmer == canonical(kmer, k);
}

} // namespace flat_index

#endif // KMER_HPP
