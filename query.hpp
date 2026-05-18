#ifndef QUERY_HPP
#define QUERY_HPP

#include <vector>
#include <cstdint>
#include "flat_assembler.hpp"
#include "mphf.hpp"

namespace flat_index {

struct RefHit {
    uint16_t ref_id;
    uint32_t ref_pos;
    bool     orient; // false: +, true: -
};

struct FullIndex {
    FlatIndex   flat_idx;
    StaticMPHF  mphf;
    std::vector<uint64_t> pos_table; // slot -> (block_byte_off << 8 | local_kmer_in_block)
};

std::vector<RefHit> query(const FullIndex& idx, const std::string& kmer);

struct Alignment {
    uint16_t unitig_id;
    bool     orient;
    uint32_t start_kmer_idx;
    uint32_t len_kmers;
};

std::vector<Alignment> align_read(const FullIndex& idx, const std::string& read);

FullIndex build_full_index(const std::vector<InUnitig>& unitigs,
                          const std::vector<InOccurrence>& occs,
                          const std::vector<InEdge>& edges,
                          int k);

} // namespace flat_index

#endif // QUERY_HPP
