#include <iostream>
#include <vector>
#include <string>
#include <cassert>
#include "data_model.hpp"
#include "query.hpp"

using namespace flat_index;

void test_tiny_fixture() {
    int k = 7;
    
    // 3 unitigs
    std::vector<InUnitig> unitigs = {
        {0, "ACGTACGT"}, // kmers: ACGTACG, CGTACGT
        {1, "TACGTAC"}, // kmers: TACGTAC
        {2, "GGGGGGG"}  // kmers: GGGGGGG
    };
    
    // Occurrences
    std::vector<InOccurrence> occs = {
        {0, 0, 100, 0, 2, '+'}, // unitig 0 at ref 0 pos 100
        {1, 0, 107, 0, 1, '+'}, // unitig 1 at ref 0 pos 107
        {2, 0, 200, 0, 1, '-'}  // unitig 2 at ref 0 pos 200, RC
    };
    
    // Edges
    std::vector<InEdge> edges = {
        {0, '+', 1, '+', 7} // unitig 0 -> unitig 1, overlap 7 (CGTACGT)
    };
    
    FullIndex idx = build_full_index(unitigs, occs, edges, k);
    
    std::cout << "Flat array size: " << idx.flat_idx.flat.size() << " bytes" << std::endl;
    
    // Query ACGTACG (unitig 0, offset 0)
    auto hits1 = query(idx, "ACGTACG");
    std::cout << "Hits for ACGTACG: " << hits1.size() << std::endl;
    for (const auto& h : hits1) {
        std::cout << "  ref=" << h.ref_id << " pos=" << h.ref_pos << " orient=" << (h.orient ? '-' : '+') << std::endl;
    }
    // ACGTACG is in unitig 0 at offset 0 (kc=ACGTACG)
    // CGTACGT is in unitig 0 at offset 1 (kc=ACGTACG)
    // Both map to kc=ACGTACG. Querying "ACGTACG" should find BOTH.
    assert(hits1.size() == 2);
    
    // Query CGTACGT (unitig 0, offset 1)
    auto hits2 = query(idx, "CGTACGT");
    std::cout << "Hits for CGTACGT: " << hits2.size() << std::endl;
    for (const auto& h : hits2) {
        std::cout << "  ref=" << h.ref_id << " pos=" << h.ref_pos << " orient=" << (h.orient ? '-' : '+') << std::endl;
    }
    assert(hits2.size() == 2);
    
    // Query TACGTAC (unitig 1, offset 0)
    auto hits3 = query(idx, "TACGTAC");
    assert(hits3.size() == 1);
    assert(hits3[0].ref_pos == 107);
    
    // Query GGGGGGG (unitig 2, offset 0, RC)
    // GGGGGGG RC is CCCCCCC
    auto hits4 = query(idx, "GGGGGGG");
    std::cout << "Hits for GGGGGGG: " << hits4.size() << std::endl;
    assert(hits4.size() == 1);
    assert(hits4[0].ref_pos == 200);
    assert(hits4[0].orient == true); // '-'
    
    // Read extension test: "ACGTACGTAC"
    // Unitig 0: ACGTACGT
    // Unitig 1: TACGTAC
    // Overlap: CGTACG (length 6)
    // Read: ACGTACGTAC
    // Path: U0 + 'A' (from U1[6]) + 'C' (from U1[7])? Wait.
    // U0: ACGTACGT
    // U1: TACGTAC
    // Read: ACGTACGTAC
    // bases: ACGTACGT + AC
    // Extension base should be 'A' (U1[6])
    
    auto alns = align_read(idx, "ACGTACGTAC");
    std::cout << "Alignments for ACGTACGTAC: " << alns.size() << std::endl;
    for (const auto& a : alns) {
        std::cout << "  unitig=" << a.unitig_id << " len=" << a.len_kmers << " orient=" << (a.orient ? '-' : '+') << std::endl;
    }
    // Should be U0 (length 2 kmers) then jump to U1 (length 2 kmers)
    // ACGTACG, CGTACGT (U0)
    // TACGTAC (U1)
    // ACGTACGTAC (len 10) -> 4 kmers: ACGTACG, CGTACGT, GTACGTA, TACGTAC?
    // Wait, U0+U1 (overlap 6): ACGTACGT + TACGTAC -> ACGTACGTAC.
    // kmers of ACGTACGTAC: ACGTACG, CGTACGT, GTACGTA, TACGTAC.
    assert(alns.size() == 2);
    assert(alns[0].unitig_id == 0);
    assert(alns[0].len_kmers == 2);
    assert(alns[1].unitig_id == 1);
    // Unitig 1 is TACGTAC, which is 7 bases. k=7. So it has ONLY 1 k-mer.
    assert(alns[1].len_kmers == 1);
    
    std::cout << "All tiny fixture tests passed!" << std::endl;
}

int main() {
    try {
        test_tiny_fixture();
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}
