// build_unitigs.cpp
//
// Self-contained compacted de Bruijn graph builder for the flat-index
// prototype. Reads a FASTA file, extracts canonical k-mers (k <= 31), traces
// maximal non-branching paths in the canonical dBG, and writes:
//
//   <out_prefix>.unitigs.tsv : id<TAB>length<TAB>sequence
//   <out_prefix>.edges.tsv   : from_id<TAB>from_orient<TAB>to_id<TAB>to_orient<TAB>(k-1)M
//
// Algorithm summary (see README "Builder algorithm" section for the full
// derivation):
//   - 2-bit base encoding A=0 C=1 G=2 T=3, each k-mer fits in a uint64_t.
//   - "Node" = canonical k-mer u = min(u, RC(u)). Each canonical k-mer has
//     two orientations in the directed graph: (u,+) reads u as written,
//     (u,-) reads RC(u).
//   - Forward neighbors of (u,o): up to 4, by appending each of A/C/G/T to
//     the end of (o == '+' ? u : RC(u)) and looking up the canonical form.
//   - In-degree of (u,o) = out-degree of (u,~o) by symmetry; computed via
//     forward_neighbors only.
//   - A unitig walk-start (u,o) is one whose in-degree != 1, OR whose unique
//     predecessor has out-degree != 1.
//   - For each unvisited canonical k-mer, try (u,+) then (u,-) as a walk
//     start; trace forward while current node has out-degree 1 AND the
//     successor has in-degree 1; mark every canonical k-mer in the walk as
//     visited so the RC-twin walk is skipped.
//   - Anything still unvisited after the linear pass lies on a pure cycle;
//     trace once with walk_cycle.
//   - For each unitig V we register two entry signatures:
//       (V.start_kmer, V.start_orient)              -> enter V at '+'
//       (V.end_kmer,   opposite(V.end_orient))      -> enter V at '-'
//     Edges are emitted from both ends: forward of V's emit-end (from='+')
//     and forward of V's emit-start opposite-orient (from='-').
//
// Memory: ~32-48 bytes per unique canonical k-mer in the std::unordered_set.
// Comfortably handles inputs up to ~100 Mb of unique k-mers on a laptop;
// scale up the memory if you point this at a full mammalian genome.
//
// Build:  make            (top-level Makefile produces bin/build_unitigs)
// Usage:  bin/build_unitigs <input.fasta> <output_prefix> [k]

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace {

constexpr int kMaxK = 31;

inline int base_to_2bit(char c) {
    switch (c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default:            return -1;
    }
}

inline char bit2_to_base(int b) {
    static constexpr char kAlphabet[4] = {'A', 'C', 'G', 'T'};
    return kAlphabet[b & 3];
}

inline std::uint64_t reverse_complement(std::uint64_t kmer, int k) {
    std::uint64_t rc = 0;
    for (int i = 0; i < k; ++i) {
        rc = (rc << 2) | ((kmer & 0x3ULL) ^ 0x3ULL);
        kmer >>= 2;
    }
    return rc;
}

inline std::string kmer_to_string(std::uint64_t kmer, int k) {
    std::string s(static_cast<std::size_t>(k), 'A');
    for (int i = k - 1; i >= 0; --i) {
        s[static_cast<std::size_t>(i)] = bit2_to_base(static_cast<int>(kmer & 0x3ULL));
        kmer >>= 2;
    }
    return s;
}

const auto opposite_orient = [](char o) -> char {
    return o == '+' ? '-' : '+';
};

struct PairHash {
    std::size_t operator()(const std::pair<std::uint64_t, char>& p) const noexcept {
        const std::uint64_t a = p.first;
        const std::uint64_t b = static_cast<std::uint64_t>(static_cast<unsigned char>(p.second));
        std::uint64_t h = a * 0x9E3779B97F4A7C15ULL;
        h ^= b + 0x9E3779B97F4A7C15ULL + (h << 6) + (h >> 2);
        return static_cast<std::size_t>(h);
    }
};

struct Neighbor {
    std::uint64_t kmer_canon;
    char orient;
    int base;
};

inline std::vector<Neighbor> forward_neighbors(std::uint64_t kmer_canon,
                                               char orient,
                                               int k,
                                               std::uint64_t mask,
                                               const std::unordered_set<std::uint64_t>& kmers) {
    std::vector<Neighbor> result;
    result.reserve(4);
    const std::uint64_t k_fwd = (orient == '+') ? kmer_canon : reverse_complement(kmer_canon, k);
    for (int b = 0; b < 4; ++b) {
        const std::uint64_t cand = ((k_fwd << 2) | static_cast<std::uint64_t>(b)) & mask;
        const std::uint64_t can = std::min(cand, reverse_complement(cand, k));
        if (kmers.find(can) != kmers.end()) {
            result.push_back({can, (cand == can) ? '+' : '-', b});
        }
    }
    return result;
}

inline bool is_walk_start(std::uint64_t u,
                          char o,
                          int k,
                          std::uint64_t mask,
                          const std::unordered_set<std::uint64_t>& kmers) {
    const auto preds = forward_neighbors(u, opposite_orient(o), k, mask, kmers);
    if (preds.size() != 1) {
        return true;
    }
    if (preds[0].kmer_canon == u) {
        // The unique predecessor is u itself in opposite orientation: this is
        // a canonical self-loop edge (X, -) -> (X, +). In canonical compacted
        // dBG semantics u is its own unitig, so treat this as a walk-start.
        return true;
    }
    const char pred_orient = opposite_orient(preds[0].orient);
    const auto pred_succs = forward_neighbors(preds[0].kmer_canon, pred_orient, k, mask, kmers);
    return pred_succs.size() != 1;
}

struct Unitig {
    std::string seq;
    std::uint64_t start_kmer = 0;
    char start_orient = '+';
    std::uint64_t end_kmer = 0;
    char end_orient = '+';
};

struct WalkResult {
    Unitig unitig;
    std::vector<std::uint64_t> visited_canons;
};

inline WalkResult walk_from(std::uint64_t u,
                            char o,
                            int k,
                            std::uint64_t mask,
                            const std::unordered_set<std::uint64_t>& kmers) {
    WalkResult r;
    r.unitig.start_kmer = u;
    r.unitig.start_orient = o;
    const std::uint64_t k_fwd = (o == '+') ? u : reverse_complement(u, k);
    r.unitig.seq = kmer_to_string(k_fwd, k);
    r.visited_canons.push_back(u);

    // Per-walk set of canonicals already in this unitig. Required to enforce
    // canonical-compacted-dBG semantics: each canonical k-mer must appear in
    // at most one unitig, so any walk step that would revisit a canonical
    // (immediate self-loop OR a longer "U-turn" through palindromic-like
    // patterns) is a unitig boundary.
    std::unordered_set<std::uint64_t> walk_canons;
    walk_canons.insert(u);

    std::uint64_t cur = u;
    char cur_o = o;
    while (true) {
        const auto succs = forward_neighbors(cur, cur_o, k, mask, kmers);
        if (succs.size() != 1) {
            break;
        }
        const auto& nb = succs[0];
        if (walk_canons.find(nb.kmer_canon) != walk_canons.end()) {
            break;
        }
        const auto nb_preds = forward_neighbors(nb.kmer_canon, opposite_orient(nb.orient), k, mask, kmers);
        if (nb_preds.size() != 1) {
            break;
        }
        if (nb_preds[0].kmer_canon == nb.kmer_canon) {
            // The successor's unique predecessor is itself: it sits at a
            // canonical self-loop and must be its own unitig.
            break;
        }
        r.unitig.seq.push_back(bit2_to_base(nb.base));
        r.visited_canons.push_back(nb.kmer_canon);
        walk_canons.insert(nb.kmer_canon);
        cur = nb.kmer_canon;
        cur_o = nb.orient;
    }
    r.unitig.end_kmer = cur;
    r.unitig.end_orient = cur_o;
    return r;
}

inline WalkResult walk_cycle(std::uint64_t u,
                             char o,
                             int k,
                             std::uint64_t mask,
                             const std::unordered_set<std::uint64_t>& kmers) {
    WalkResult r;
    r.unitig.start_kmer = u;
    r.unitig.start_orient = o;
    const std::uint64_t k_fwd = (o == '+') ? u : reverse_complement(u, k);
    r.unitig.seq = kmer_to_string(k_fwd, k);
    r.visited_canons.push_back(u);

    std::unordered_set<std::uint64_t> walk_canons;
    walk_canons.insert(u);

    std::uint64_t cur = u;
    char cur_o = o;
    while (true) {
        const auto succs = forward_neighbors(cur, cur_o, k, mask, kmers);
        if (succs.size() != 1) {
            break;
        }
        const auto& nb = succs[0];
        if (walk_canons.find(nb.kmer_canon) != walk_canons.end()) {
            break;
        }
        r.unitig.seq.push_back(bit2_to_base(nb.base));
        r.visited_canons.push_back(nb.kmer_canon);
        walk_canons.insert(nb.kmer_canon);
        cur = nb.kmer_canon;
        cur_o = nb.orient;
    }
    r.unitig.end_kmer = cur;
    r.unitig.end_orient = cur_o;
    return r;
}

inline std::vector<std::string> read_fasta(const std::string& path) {
    std::ifstream in(path);
    if (!in.is_open()) {
        throw std::runtime_error("cannot open input FASTA: " + path);
    }
    std::vector<std::string> seqs;
    std::string line;
    std::string current;
    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }
        if (!line.empty() && line[0] == '>') {
            if (!current.empty()) {
                seqs.push_back(std::move(current));
                current.clear();
            }
            continue;
        }
        current += line;
    }
    if (!current.empty()) {
        seqs.push_back(std::move(current));
    }
    return seqs;
}

} // namespace

int main(int argc, char** argv) {
    if (argc < 3 || argc > 4) {
        std::fprintf(stderr,
                     "Usage: %s <input.fasta> <output_prefix> [k]\n"
                     "  k defaults to 31; max %d.\n",
                     argv[0], kMaxK);
        return 1;
    }
    const std::string input_path = argv[1];
    const std::string out_prefix = argv[2];
    const int k = (argc == 4) ? std::atoi(argv[3]) : 31;
    if (k < 2 || k > kMaxK) {
        std::fprintf(stderr, "error: k must be in [2, %d]\n", kMaxK);
        return 1;
    }
    const std::uint64_t mask = (k == 32)
        ? ~static_cast<std::uint64_t>(0)
        : ((static_cast<std::uint64_t>(1) << (2 * k)) - 1);

    std::vector<std::string> sequences;
    try {
        sequences = read_fasta(input_path);
    } catch (const std::exception& e) {
        std::fprintf(stderr, "[build_unitigs] %s\n", e.what());
        return 1;
    }
    std::fprintf(stderr, "[build_unitigs] read %zu sequence(s) from %s\n",
                 sequences.size(), input_path.c_str());

    std::unordered_set<std::uint64_t> kmers;
    {
        std::size_t windows_seen = 0;
        for (const auto& seq : sequences) {
            std::uint64_t fwd = 0;
            std::uint64_t rev = 0;
            int filled = 0;
            for (char c : seq) {
                const int b = base_to_2bit(c);
                if (b < 0) {
                    filled = 0;
                    fwd = 0;
                    rev = 0;
                    continue;
                }
                fwd = ((fwd << 2) | static_cast<std::uint64_t>(b)) & mask;
                rev = (rev >> 2) |
                      (static_cast<std::uint64_t>(b ^ 3) << (2 * (k - 1)));
                if (++filled >= k) {
                    kmers.insert(std::min(fwd, rev));
                    ++windows_seen;
                }
            }
        }
        std::fprintf(stderr,
                     "[build_unitigs] saw %zu kmer windows; %zu unique canonical k-mers\n",
                     windows_seen, kmers.size());
    }

    std::vector<Unitig> unitigs;
    std::unordered_set<std::uint64_t> visited;
    visited.reserve(kmers.size() * 2);

    for (std::uint64_t u : kmers) {
        if (visited.find(u) != visited.end()) {
            continue;
        }
        WalkResult r;
        if (is_walk_start(u, '+', k, mask, kmers)) {
            r = walk_from(u, '+', k, mask, kmers);
        } else if (is_walk_start(u, '-', k, mask, kmers)) {
            r = walk_from(u, '-', k, mask, kmers);
        } else {
            continue;
        }
        for (std::uint64_t c : r.visited_canons) {
            visited.insert(c);
        }
        unitigs.push_back(std::move(r.unitig));
    }

    for (std::uint64_t u : kmers) {
        if (visited.find(u) != visited.end()) {
            continue;
        }
        WalkResult r = walk_cycle(u, '+', k, mask, kmers);
        for (std::uint64_t c : r.visited_canons) {
            visited.insert(c);
        }
        unitigs.push_back(std::move(r.unitig));
    }

    std::fprintf(stderr, "[build_unitigs] built %zu unitigs\n", unitigs.size());

    using SigKey = std::pair<std::uint64_t, char>;
    std::unordered_map<SigKey, std::pair<int, char>, PairHash> entry_lookup;
    entry_lookup.reserve(unitigs.size() * 2);
    for (std::size_t i = 0; i < unitigs.size(); ++i) {
        const auto& V = unitigs[i];
        entry_lookup[{V.start_kmer, V.start_orient}] = {static_cast<int>(i), '+'};
        entry_lookup[{V.end_kmer, opposite_orient(V.end_orient)}] = {static_cast<int>(i), '-'};
    }

    struct Edge {
        int from_id;
        char from_orient;
        int to_id;
        char to_orient;
        int overlap;
    };
    std::vector<Edge> edges;

    for (std::size_t i = 0; i < unitigs.size(); ++i) {
        const auto& V = unitigs[i];
        const auto succs_fwd = forward_neighbors(V.end_kmer, V.end_orient, k, mask, kmers);
        for (const auto& s : succs_fwd) {
            const auto it = entry_lookup.find({s.kmer_canon, s.orient});
            if (it == entry_lookup.end()) {
                continue;
            }
            edges.push_back({static_cast<int>(i), '+', it->second.first, it->second.second, k - 1});
        }
        const auto succs_bwd = forward_neighbors(V.start_kmer, opposite_orient(V.start_orient), k, mask, kmers);
        for (const auto& s : succs_bwd) {
            const auto it = entry_lookup.find({s.kmer_canon, s.orient});
            if (it == entry_lookup.end()) {
                continue;
            }
            edges.push_back({static_cast<int>(i), '-', it->second.first, it->second.second, k - 1});
        }
    }

    std::fprintf(stderr, "[build_unitigs] built %zu edges\n", edges.size());

    const std::string unitigs_path = out_prefix + ".unitigs.tsv";
    const std::string edges_path = out_prefix + ".edges.tsv";

    {
        const std::filesystem::path out_path(out_prefix);
        const auto parent = out_path.parent_path();
        if (!parent.empty()) {
            std::error_code ec;
            std::filesystem::create_directories(parent, ec);
            if (ec) {
                std::fprintf(stderr,
                             "[build_unitigs] error: cannot create output directory %s: %s\n",
                             parent.string().c_str(), ec.message().c_str());
                return 1;
            }
        }
    }

    std::ofstream u_out(unitigs_path);
    std::ofstream e_out(edges_path);
    if (!u_out.is_open() || !e_out.is_open()) {
        std::fprintf(stderr, "[build_unitigs] error: cannot write outputs under %s\n",
                     out_prefix.c_str());
        return 1;
    }

    u_out << "id\tlength\tsequence\n";
    for (std::size_t i = 0; i < unitigs.size(); ++i) {
        u_out << i << '\t' << unitigs[i].seq.size() << '\t' << unitigs[i].seq << '\n';
    }

    e_out << "from_id\tfrom_orient\tto_id\tto_orient\toverlap\n";
    for (const auto& e : edges) {
        e_out << e.from_id << '\t' << e.from_orient << '\t'
              << e.to_id << '\t' << e.to_orient << '\t' << e.overlap << "M\n";
    }

    std::size_t total_bp = 0;
    std::size_t min_len = std::numeric_limits<std::size_t>::max();
    std::size_t max_len = 0;
    for (const auto& V : unitigs) {
        total_bp += V.seq.size();
        if (V.seq.size() < min_len) {
            min_len = V.seq.size();
        }
        if (V.seq.size() > max_len) {
            max_len = V.seq.size();
        }
    }
    if (unitigs.empty()) {
        min_len = 0;
    }
    const double avg_len = unitigs.empty()
        ? 0.0
        : static_cast<double>(total_bp) / static_cast<double>(unitigs.size());

    std::fprintf(stderr,
                 "[build_unitigs] unitigs=%zu edges=%zu total_bp=%zu min=%zu max=%zu avg=%.2f\n"
                 "[build_unitigs] wrote %s\n"
                 "[build_unitigs] wrote %s\n",
                 unitigs.size(), edges.size(), total_bp, min_len, max_len, avg_len,
                 unitigs_path.c_str(), edges_path.c_str());

    return 0;
}
