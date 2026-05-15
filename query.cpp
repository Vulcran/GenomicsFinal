#include "query.hpp"
#include "kmer.hpp"
#include "dna_pack.hpp"
#include "flat_block.hpp"
#include "metrics.hpp"
#include <unordered_map>
#include <iostream>

namespace flat_index {

FullIndex build_full_index(const std::vector<InUnitig>& unitigs,
                          const std::vector<InOccurrence>& occs,
                          const std::vector<InEdge>& edges,
                          int k) {
    FullIndex idx;
    idx.flat_idx = assemble(unitigs, occs, edges, k);
    
    std::unordered_map<uint64_t, uint64_t> kmer_to_pos;
    kmer_to_pos.reserve((idx.flat_idx.flat.size() / BLOCK_SIZE) * 128);
    
    size_t current_offset = 0;
    size_t head_offset = 0;
    uint16_t last_uid = 0xFFFF;
    size_t total_blocks = idx.flat_idx.flat.size() / BLOCK_SIZE;
    size_t blocks_processed = 0;

    while (current_offset < idx.flat_idx.flat.size()) {
        if (blocks_processed % 1000 == 0 || blocks_processed == total_blocks - 1) {
            std::cout << "\rCollecting k-mers: " << (blocks_processed + 1) << " / " << total_blocks << std::flush;
        }
        const std::byte* block_ptr = &idx.flat_idx.flat[current_offset];
        BlockHeader h = BlockHeader::unpack(block_ptr);
        
        if (h.unitig_id != last_uid) {
            head_offset = current_offset;
            last_uid = h.unitig_id;
        }
        
        const std::byte* dna_ptr = block_ptr + DNA_REGION_START;
        uint64_t k_enc = 0;
        uint64_t mask = (k == 32) ? ~0ULL : (1ULL << (2 * k)) - 1;

        if (h.dna_len_kmers > 0) {
            k_enc = read_kmer(dna_ptr, 0, k);
        }

        for (uint8_t i = 0; i < h.dna_len_kmers; ++i) {
            if (i > 0) {
                size_t next_bit_off = (i + k - 1) * 2;
                uint8_t b = static_cast<uint8_t>(dna_ptr[next_bit_off / 8]);
                uint8_t next_base = (b >> (6 - (next_bit_off % 8))) & 0x03;
                k_enc = ((k_enc << 2) | next_base) & mask;
            }
            uint64_t kc = canonical(k_enc, k);
            
            // Optimized Pufferfish packing: [0-31] head_block_off, [32-63] kmer_idx_in_unitig
            uint32_t kmer_idx_in_unitig = (current_offset - head_offset) / BLOCK_SIZE * 128 + i;
            uint64_t pos_entry = (static_cast<uint64_t>(kmer_idx_in_unitig) << 32) | static_cast<uint32_t>(head_offset);
            kmer_to_pos[kc] = pos_entry;
        }
        
        current_offset += BLOCK_SIZE;
        blocks_processed++;
    }
    std::cout << std::endl;
    
    std::vector<uint64_t> keys;
    keys.reserve(kmer_to_pos.size());
    for (const auto& kv : kmer_to_pos) {
        keys.push_back(kv.first);
    }
    
    idx.mphf.build(keys);
    idx.pos_table.resize(idx.mphf.table_size());
    
    size_t keys_processed = 0;
    size_t total_keys = kmer_to_pos.size();
    for (const auto& kv : kmer_to_pos) {
        if (keys_processed % 10000 == 0 || keys_processed == total_keys - 1) {
            std::cout << "\rAssigning positions: " << (keys_processed + 1) << " / " << total_keys << std::flush;
        }
        size_t slot = idx.mphf.lookup(kv.first);
        if (slot >= idx.pos_table.size()) idx.pos_table.resize(slot + 1);
        idx.pos_table[slot] = kv.second;
        keys_processed++;
    }
    std::cout << std::endl;
    
    return idx;
}

std::vector<RefHit> query(const FullIndex& idx, const std::string& kmer_str) {
    int k = idx.flat_idx.k;
    if (kmer_str.length() != (size_t)k) return {};
    
    uint64_t k_enc = encode(kmer_str);
    uint64_t kc = canonical(k_enc, k);
    
    size_t slot = idx.mphf.lookup(kc);
    if (slot >= idx.pos_table.size()) return {};
    
    g_cache_counter.record_miss(); // Access pos_table
    uint64_t pos = idx.pos_table[slot];
    uint32_t head_off = static_cast<uint32_t>(pos & 0xFFFFFFFFULL);
    uint32_t kmer_idx_in_unitig = static_cast<uint32_t>(pos >> 32);
    
    uint32_t block_off = head_off + (kmer_idx_in_unitig / 128) * BLOCK_SIZE;
    uint8_t loc_kmers = static_cast<uint8_t>(kmer_idx_in_unitig % 128);
    
    if (block_off >= idx.flat_idx.flat.size()) return {};
    g_cache_counter.record_miss(); // Access flat array block
    const std::byte* block_ptr = &idx.flat_idx.flat[block_off];
    uint64_t k_read = read_kmer(block_ptr + DNA_REGION_START, loc_kmers * 2, k);
    if (canonical(k_read, k) != kc) return {};
    
    bool query_is_rc = (k_enc != k_read);

    const std::byte* head_ptr = &idx.flat_idx.flat[head_off];
    if (head_off != block_off) {
        g_cache_counter.record_miss(); // Access head block (different from k-mer block)
    }
    BlockHeader h = BlockHeader::unpack(head_ptr);
    
    size_t utab_start = DNA_REGION_START + DNA_REGION_SIZE;
    std::vector<RefHit> all_hits;
    for (int i = 0; i < h.n_utab; ++i) {
        UTabRow row = UTabRow::unpack(head_ptr + utab_start + i * 8);
        RefHit hit;
        hit.ref_id = row.ref_id;
        hit.orient = (row.orient != query_is_rc);
        if (!hit.orient) hit.ref_pos = row.ref_pos + (kmer_idx_in_unitig - row.entry_off);
        else hit.ref_pos = row.ref_pos + (row.entry_off - kmer_idx_in_unitig);
        all_hits.push_back(hit);
    }
    return all_hits;
}

std::vector<Alignment> align_read(const FullIndex& idx, const std::string& read) {
    int k = idx.flat_idx.k;
    if (read.length() < (size_t)k) return {};
    
    std::string first_kmer = read.substr(0, k);
    uint64_t k_enc = encode(first_kmer);
    uint64_t kc = canonical(k_enc, k);
    
    size_t slot = idx.mphf.lookup(kc);
    if (slot >= idx.pos_table.size()) return {};
    
    uint64_t pos = idx.pos_table[slot];
    uint32_t block_off = static_cast<uint32_t>(pos & 0x00FFFFFFFFFFFFFFULL);
    uint8_t loc_kmers = static_cast<uint8_t>(pos >> 56);
    
    if (block_off >= idx.flat_idx.flat.size()) return {};
    const std::byte* block_ptr = &idx.flat_idx.flat[block_off];
    uint64_t k_read = read_kmer(block_ptr + DNA_REGION_START, loc_kmers * 2, k);
    if (canonical(k_read, k) != kc) return {};
    
    bool current_orient = (k_enc != k_read);
    uint32_t current_block_off = block_off;
    uint8_t current_loc_kmers = loc_kmers;
    
    BlockHeader h = BlockHeader::unpack(block_ptr);
    Alignment curr_aln;
    curr_aln.unitig_id = h.unitig_id;
    curr_aln.orient = current_orient;
    
    uint32_t head_off = current_block_off;
    while (head_off >= BLOCK_SIZE) {
        BlockHeader prev_h = BlockHeader::unpack(&idx.flat_idx.flat[head_off - BLOCK_SIZE]);
        if (prev_h.unitig_id != h.unitig_id || (prev_h.flags & 0x01) == 0) break;
        head_off -= BLOCK_SIZE;
    }
    uint32_t kmers_per_block = 128;
    curr_aln.start_kmer_idx = (current_block_off - head_off) / BLOCK_SIZE * kmers_per_block + current_loc_kmers;
    curr_aln.len_kmers = 1;
    
    std::vector<Alignment> alns;
    for (size_t i = 1; i + k <= read.length(); ++i) {
        char next_base = read[i + k - 1];
        bool extended = false;
        
        if (!current_orient) { // Forward
            if (current_loc_kmers + 1 < h.dna_len_kmers) {
                current_loc_kmers++;
                extended = true;
            } else if (h.flags & 0x01) {
                current_block_off += BLOCK_SIZE;
                current_loc_kmers = 0;
                block_ptr = &idx.flat_idx.flat[current_block_off];
                h = BlockHeader::unpack(block_ptr);
                extended = true;
            } else {
                size_t etab_start = DNA_REGION_START + DNA_REGION_SIZE + h.n_utab * 8;
                for (int e = 0; e < h.n_etab; ++e) {
                    ETabRow er = ETabRow::unpack(block_ptr + etab_start + e * 4);
                    if (bits_to_base(er.ext_base) == next_base && er.from_orient == false) {
                        alns.push_back(curr_aln);
                        current_block_off = er.to_block_offset_div_64 * 64;
                        current_loc_kmers = 0;
                        current_orient = er.to_orient;
                        block_ptr = &idx.flat_idx.flat[current_block_off];
                        h = BlockHeader::unpack(block_ptr);
                        curr_aln.unitig_id = h.unitig_id;
                        curr_aln.orient = current_orient;
                        curr_aln.start_kmer_idx = 0;
                        curr_aln.len_kmers = 1;
                        extended = true;
                        break;
                    }
                }
            }
        } else { // RC
            if (current_loc_kmers > 0) {
                current_loc_kmers--;
                extended = true;
            } else {
                bool has_prev = false;
                if (current_block_off >= BLOCK_SIZE) {
                    BlockHeader prev_h = BlockHeader::unpack(&idx.flat_idx.flat[current_block_off - BLOCK_SIZE]);
                    if (prev_h.unitig_id == h.unitig_id && (prev_h.flags & 0x01)) {
                        current_block_off -= BLOCK_SIZE;
                        block_ptr = &idx.flat_idx.flat[current_block_off];
                        h = BlockHeader::unpack(block_ptr);
                        current_loc_kmers = h.dna_len_kmers - 1;
                        extended = true;
                        has_prev = true;
                    }
                }
                if (!has_prev) {
                    size_t etab_start = DNA_REGION_START + DNA_REGION_SIZE + h.n_utab * 8;
                    for (int e = 0; e < h.n_etab; ++e) {
                        ETabRow er = ETabRow::unpack(block_ptr + etab_start + e * 4);
                        if (bits_to_base(er.ext_base) == next_base && er.from_orient == true) {
                            alns.push_back(curr_aln);
                            uint32_t target_head_off = er.to_block_offset_div_64 * 64;
                            current_orient = er.to_orient;
                            if (current_orient) {
                                uint32_t t_off = target_head_off;
                                BlockHeader th = BlockHeader::unpack(&idx.flat_idx.flat[t_off]);
                                while (th.flags & 0x01) {
                                    t_off += BLOCK_SIZE;
                                    th = BlockHeader::unpack(&idx.flat_idx.flat[t_off]);
                                }
                                current_block_off = t_off;
                                block_ptr = &idx.flat_idx.flat[current_block_off];
                                h = th;
                                current_loc_kmers = h.dna_len_kmers - 1;
                            } else {
                                current_block_off = target_head_off;
                                block_ptr = &idx.flat_idx.flat[current_block_off];
                                h = BlockHeader::unpack(block_ptr);
                                current_loc_kmers = 0;
                            }
                            curr_aln.unitig_id = h.unitig_id;
                            curr_aln.orient = current_orient;
                            if (current_orient) {
                                uint32_t t_head = target_head_off;
                                uint32_t total_kmers = 0;
                                uint32_t t_curr = t_head;
                                while(true) {
                                    BlockHeader th = BlockHeader::unpack(&idx.flat_idx.flat[t_curr]);
                                    total_kmers += th.dna_len_kmers;
                                    if (!(th.flags & 0x01)) break;
                                    t_curr += BLOCK_SIZE;
                                }
                                curr_aln.start_kmer_idx = total_kmers - 1;
                            } else {
                                curr_aln.start_kmer_idx = 0;
                            }
                            curr_aln.len_kmers = 1;
                            extended = true;
                            break;
                        }
                    }
                }
            }
        }
        
        if (extended) {
            uint64_t k_read_ext = read_kmer(block_ptr + DNA_REGION_START, current_loc_kmers * 2, k);
            uint64_t k_expected = encode(read.substr(i, k));
            if (current_orient) k_read_ext = reverse_complement(k_read_ext, k);
            if (k_read_ext == k_expected) curr_aln.len_kmers++;
            else break;
        } else break;
    }
    alns.push_back(curr_aln);
    return alns;
}

} // namespace flat_index
