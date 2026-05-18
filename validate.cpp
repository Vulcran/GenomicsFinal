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
    // The POS table stores ONE entry per canonical k-mer. ACGTACG and CGTACGT share
    // canonical ACGTACG; the build loop stores the last-seen position (local kmer 1 =
    // CGTACGT). Querying "ACGTACG" finds that position, so 1 hit.
    auto hits1 = query(idx, "ACGTACG");
    std::cout << "Hits for ACGTACG: " << hits1.size() << std::endl;
    for (const auto& h : hits1) {
        std::cout << "  ref=" << h.ref_id << " pos=" << h.ref_pos << " orient=" << (h.orient ? '-' : '+') << std::endl;
    }
    assert(hits1.size() == 1);
    assert(hits1[0].ref_id == 0);
    assert(hits1[0].ref_pos == 101); // stored CGTACGT is at pos 101; query is its RC → orient '-'
    assert(hits1[0].orient == true); // '-'

    // Query CGTACGT — same canonical, same stored slot, forward match this time
    auto hits2 = query(idx, "CGTACGT");
    std::cout << "Hits for CGTACGT: " << hits2.size() << std::endl;
    for (const auto& h : hits2) {
        std::cout << "  ref=" << h.ref_id << " pos=" << h.ref_pos << " orient=" << (h.orient ? '-' : '+') << std::endl;
    }
    assert(hits2.size() == 1);
    assert(hits2[0].ref_pos == 101);
    assert(hits2[0].orient == false); // '+'

    // Query TACGTAC (unitig 1, offset 0)
    auto hits3 = query(idx, "TACGTAC");
    assert(hits3.size() == 1);
    assert(hits3[0].ref_pos == 107);

    // Query GGGGGGG (unitig 2, offset 0, unitig placed in '-' orientation in ref)
    auto hits4 = query(idx, "GGGGGGG");
    std::cout << "Hits for GGGGGGG: " << hits4.size() << std::endl;
    assert(hits4.size() == 1);
    assert(hits4[0].ref_pos == 200);
    assert(hits4[0].orient == true); // '-'

    // Read extension: "ACGTACGTAC" seeds on ACGTACG.
    // The POS entry points to local kmer 1 (CGTACGT), so align_read starts in RC
    // orientation and walks backward through unitig 0. The U0→U1 edge has from_orient='+'
    // and RC mode looks for from_orient='-', so the walk stops at unitig 0.
    auto alns = align_read(idx, "ACGTACGTAC");
    std::cout << "Alignments for ACGTACGTAC: " << alns.size() << std::endl;
    for (const auto& a : alns) {
        std::cout << "  unitig=" << a.unitig_id << " len=" << a.len_kmers << " orient=" << (a.orient ? '-' : '+') << std::endl;
    }
    assert(alns.size() == 1);
    assert(alns[0].unitig_id == 0);
    assert(alns[0].len_kmers == 2);
    
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
