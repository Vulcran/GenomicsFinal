#include "block_builder.hpp"
#include "dna_pack.hpp"
#include <algorithm>
#include <iostream>

namespace flat_index {

PackedBlocks pack_unitig(const InUnitig& u,
                         const std::vector<InOccurrence>& utab_rows,
                         const std::vector<InEdge>& etab_rows,
                         int k) {
    PackedBlocks pb;
    
    // We target DNA_REGION_SIZE = 40 bytes for DNA to leave space for UTAB/ETAB.
    // D = 40 bytes = 320 bits = 160 bases.
    // Each block covers BASES_PER_BLOCK = 128 new bases.
    // To allow read_kmer to work without crossing blocks, we store up to k-1 additional 
    // bases from the next block(s) if they belong to a k-mer starting in this block.
    
    const size_t D = DNA_REGION_SIZE; 
    const size_t BASES_PER_BLOCK = 128; 
    
    size_t total_bases = u.seq.length();
    size_t num_blocks = (total_bases + BASES_PER_BLOCK - 1) / BASES_PER_BLOCK;
    if (num_blocks == 0) num_blocks = 1;
    
    pb.bytes.resize(num_blocks * BLOCK_SIZE, std::byte{0});
    
    for (size_t i = 0; i < num_blocks; ++i) {
        std::byte* block_ptr = &pb.bytes[i * BLOCK_SIZE];
        
        BlockHeader header;
        header.unitig_id = static_cast<uint16_t>(u.id);
        
        size_t start_base = i * BASES_PER_BLOCK;
        size_t end_base = std::min((i + 1) * BASES_PER_BLOCK, total_bases);
        
        // dna_len_kmers is how many k-mers START in this block.
        header.dna_len_kmers = (total_bases - start_base >= (size_t)k) ? 
                               (uint8_t)std::min((size_t)BASES_PER_BLOCK, total_bases - start_base - k + 1) : 0;
        
        // Special case: if this is not the last block, it might have exactly BASES_PER_BLOCK k-mers starting in it
        if (i < num_blocks - 1) {
            header.dna_len_kmers = (uint8_t)BASES_PER_BLOCK;
        }

        size_t dna_end_base = std::min(start_base + D * 4, total_bases);
        size_t block_dna_bases = dna_end_base - start_base;
        
        std::vector<InEdge> head_edges, tail_edges;
        for (const auto& et : etab_rows) {
            if (et.from_orient == '-') head_edges.push_back(et);
            else tail_edges.push_back(et);
        }

        header.n_utab = (i == 0) ? (uint8_t)utab_rows.size() : 0;
        if (num_blocks == 1) {
            header.n_etab = (uint8_t)etab_rows.size();
        } else {
            if (i == 0) header.n_etab = (uint8_t)head_edges.size();
            else if (i == num_blocks - 1) header.n_etab = (uint8_t)tail_edges.size();
            else header.n_etab = 0;
        }
        
        header.flags = (i < num_blocks - 1) ? 0x01 : 0x00; // 0x01 = has_continuation
        header.next_block_offset_div_64 = (i < num_blocks - 1) ? 1 : 0; 
        
        header.pack(block_ptr);
        
        // Pack DNA
        std::string_view subseq = std::string_view(u.seq).substr(start_base, block_dna_bases);
        pack_dna(subseq, block_ptr + DNA_REGION_START);
        
        if (i == 0 || (i == num_blocks - 1 && header.n_etab > 0)) {
            size_t current_off = DNA_REGION_START + D;
            if (i == 0) {
                // Pack UTAB
                for (const auto& ut : utab_rows) {
                    if (current_off + 8 > BLOCK_SIZE) break;
                    UTabRow row;
                    row.ref_id = ut.ref_id;
                    row.ref_pos = ut.ref_pos;
                    row.entry_off = ut.entry_off;
                    row.walk_len = ut.walk_len;
                    row.orient = (ut.orient == '-');
                    row.pack(block_ptr + current_off);
                    current_off += 8;
                }
                // Pack Head Edges
                const auto& edges_to_pack = (num_blocks == 1) ? etab_rows : head_edges;
                for (const auto& et : edges_to_pack) {
                    if (current_off + 4 > BLOCK_SIZE) break;
                    ETabRow row;
                    row.to_block_offset_div_64 = 0;
                    row.ext_base = 0;
                    row.from_orient = (et.from_orient == '-');
                    row.to_orient = (et.to_orient == '-');
                    row.pack(block_ptr + current_off);
                    current_off += 4;
                }
            } else if (i == num_blocks - 1) {
                // Pack Tail Edges (UTAB is already skipped because n_utab=0)
                for (const auto& et : tail_edges) {
                    if (current_off + 4 > BLOCK_SIZE) break;
                    ETabRow row;
                    row.to_block_offset_div_64 = 0;
                    row.ext_base = 0;
                    row.from_orient = (et.from_orient == '-');
                    row.to_orient = (et.to_orient == '-');
                    row.pack(block_ptr + current_off);
                    current_off += 4;
                }
            }
        }
    }
    
    return pb;
}

} // namespace flat_index
