#include "flat_assembler.hpp"
#include "block_builder.hpp"
#include "kmer.hpp"
#include <unordered_map>
#include <map>
#include <iostream>

namespace flat_index {

FlatIndex assemble(const std::vector<InUnitig>& unitigs,
                   const std::vector<InOccurrence>& occs,
                   const std::vector<InEdge>& edges,
                   int k) {
    FlatIndex fi;
    fi.k = k;

    std::unordered_map<uint32_t, std::vector<InOccurrence>> unitig_to_occs;
    for (const auto& occ : occs) {
        unitig_to_occs[occ.unitig_id].push_back(occ);
    }

    std::unordered_map<uint32_t, std::vector<InEdge>> unitig_to_edges;
    for (const auto& edge : edges) {
        unitig_to_edges[edge.from_id].push_back(edge);
    }

    std::unordered_map<uint32_t, uint32_t> unitig_id_to_byte_offset;
    std::vector<PackedBlocks> all_packed;

    // First pass: Pack unitigs and record offsets
    for (size_t i = 0; i < unitigs.size(); ++i) {
        if (i % 100 == 0 || i == unitigs.size() - 1) {
            std::cout << "\rPacking unitigs: " << (i + 1) << " / " << unitigs.size() << std::flush;
        }
        const auto& u = unitigs[i];
        unitig_id_to_byte_offset[u.id] = static_cast<uint32_t>(fi.flat.size());
        PackedBlocks pb = pack_unitig(u, unitig_to_occs[u.id], unitig_to_edges[u.id], k);
        
        // Fix continuation block offsets (absolute in flat array)
        size_t num_blocks = pb.bytes.size() / BLOCK_SIZE;
        for (size_t i = 0; i < num_blocks; ++i) {
            if (i < num_blocks - 1) {
                BlockHeader h = BlockHeader::unpack(&pb.bytes[i * BLOCK_SIZE]);
                h.next_block_offset_div_64 = (unitig_id_to_byte_offset[u.id] + (i + 1) * BLOCK_SIZE) / 64;
                h.pack(&pb.bytes[i * BLOCK_SIZE]);
            }
        }
        
        fi.flat.insert(fi.flat.end(), pb.bytes.begin(), pb.bytes.end());
        all_packed.push_back(std::move(pb));
    }
    std::cout << std::endl;

    // Second pass: back-patch ETAB
    std::unordered_map<uint32_t, const InUnitig*> id_to_unitig;
    for (const auto& u : unitigs) {
        id_to_unitig[u.id] = &u;
    }

    size_t current_unitig_offset = 0;
    for (size_t u_idx = 0; u_idx < unitigs.size(); ++u_idx) {
        if (u_idx % 100 == 0 || u_idx == unitigs.size() - 1) {
            std::cout << "\rPatching edges: " << (u_idx + 1) << " / " << unitigs.size() << std::flush;
        }
        const auto& u = unitigs[u_idx];
        const auto& u_edges = unitig_to_edges[u.id];
        const auto& pb = all_packed[u_idx];
        size_t num_blocks = pb.bytes.size() / BLOCK_SIZE;
        
        if (!u_edges.empty()) {
            auto patch_edges = [&](std::byte* b_ptr, const std::vector<InEdge>& edges_to_patch) {
                BlockHeader bh = BlockHeader::unpack(b_ptr);
                size_t e_start = DNA_REGION_START + DNA_REGION_SIZE + bh.n_utab * 8;
                for (size_t e_idx = 0; e_idx < edges_to_patch.size(); ++e_idx) {
                    const auto& edge = edges_to_patch[e_idx];
                    ETabRow row = ETabRow::unpack(b_ptr + e_start + e_idx * 4);
                    row.to_block_offset_div_64 = unitig_id_to_byte_offset[edge.to_id] / 64;
                    const auto& target_u = *id_to_unitig[edge.to_id];
                    if (edge.to_orient == '+') {
                        if (edge.overlap < target_u.seq.length())
                            row.ext_base = base_to_bits(target_u.seq[edge.overlap]);
                    } else {
                        if (edge.overlap < target_u.seq.length()) {
                            char b = target_u.seq[target_u.seq.length() - edge.overlap - 1];
                            row.ext_base = base_to_bits(b) ^ 0x03;
                        }
                    }
                    row.pack(b_ptr + e_start + e_idx * 4);
                }
            };

            if (num_blocks == 1) {
                patch_edges(&fi.flat[current_unitig_offset], u_edges);
            } else {
                std::vector<InEdge> head_edges, tail_edges;
                for (const auto& et : u_edges) {
                    if (et.from_orient == '-') head_edges.push_back(et);
                    else tail_edges.push_back(et);
                }
                if (!head_edges.empty()) patch_edges(&fi.flat[current_unitig_offset], head_edges);
                if (!tail_edges.empty()) patch_edges(&fi.flat[current_unitig_offset + (num_blocks - 1) * BLOCK_SIZE], tail_edges);
            }
        }
        
        current_unitig_offset += pb.bytes.size();
    }
    std::cout << std::endl;

    return fi;
}

} // namespace flat_index
