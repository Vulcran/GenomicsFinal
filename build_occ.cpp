// build_occ.cpp
//
// Builds the unitig occurrence table (UTAB / "ContigTable" in Pufferfish) for
// the flat-index prototype.
//
// Inputs:
//   <input.fasta>                the reference FASTA originally given to
//                                build_unitigs
//   <out_prefix>.unitigs.tsv     produced by build_unitigs
//
// Outputs:
//   <out_prefix>.occ.tsv   columns: unitig_id, ref_id, ref_pos, entry_off,
//                                   orient, walk_len
//   <out_prefix>.refs.tsv  columns: ref_id, name, length
//                          (kept alongside so downstream code can look up a
//                           FASTA name from the integer ref_id)
//
// Build:  make            (top-level Makefile produces bin/build_occ)
// Usage:  bin/build_occ <input.fasta> <out_prefix> [k]

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
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

inline char complement_base(char c) {
    switch (c) {
        case 'A': case 'a': return 'T';
        case 'C': case 'c': return 'G';
        case 'G': case 'g': return 'C';
        case 'T': case 't': return 'A';
        default:            return 'N';
    }
}

inline std::string reverse_complement_string(const std::string& s) {
    std::string out(s.size(), 'A');
    for (std::size_t i = 0; i < s.size(); ++i) {
        out[i] = complement_base(s[s.size() - 1 - i]);
    }
    return out;
}

struct KmerLoc {
    int unitig_id;
    int off_local;
    char o_local;
};

inline std::vector<std::string> read_fasta_with_names(const std::string& path,
                                                     std::vector<std::string>& names) {
    std::ifstream in(path);
    if (!in.is_open()) {
        throw std::runtime_error("cannot open input FASTA: " + path);
    }
    std::vector<std::string> seqs;
    std::string line;
    std::string current;
    std::string current_name;
    bool have_record = false;
    const auto flush = [&]() {
        if (have_record) {
            seqs.push_back(std::move(current));
            names.push_back(std::move(current_name));
            current.clear();
            current_name.clear();
            have_record = false;
        }
    };
    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }
        if (!line.empty() && line[0] == '>') {
            flush();
            const std::string hdr = line.substr(1);
            const std::size_t sp = hdr.find_first_of(" \t");
            current_name = (sp == std::string::npos) ? hdr : hdr.substr(0, sp);
            have_record = true;
            continue;
        }
        if (have_record) {
            current += line;
        }
    }
    flush();
    return seqs;
}

struct ParsedUnitigs {
    std::vector<std::string> seqs;
};

inline ParsedUnitigs read_unitigs_tsv(const std::string& path) {
    std::ifstream in(path);
    if (!in.is_open()) {
        throw std::runtime_error("cannot open unitigs TSV: " + path);
    }
    ParsedUnitigs out;
    std::string line;
    bool first = true;
    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }
        if (line.empty()) {
            continue;
        }
        if (first) {
            first = false;
            if (line.rfind("id\t", 0) == 0) {
                continue;
            }
        }
        std::stringstream ss(line);
        std::string id_s, len_s, seq_s;
        if (!std::getline(ss, id_s, '\t')) continue;
        if (!std::getline(ss, len_s, '\t')) continue;
        if (!std::getline(ss, seq_s, '\t')) continue;
        const int id = std::atoi(id_s.c_str());
        if (id < 0) {
            throw std::runtime_error("negative unitig id in TSV: " + path);
        }
        if (static_cast<std::size_t>(id) >= out.seqs.size()) {
            out.seqs.resize(static_cast<std::size_t>(id) + 1);
        }
        out.seqs[static_cast<std::size_t>(id)] = std::move(seq_s);
    }
    return out;
}

struct OccRow {
    int unitig_id;
    int ref_id;
    long long ref_pos;
    int entry_off;
    char orient;
    int walk_len;
};

} // namespace

int main(int argc, char** argv) {
    if (argc < 3 || argc > 4) {
        std::fprintf(stderr,
                     "Usage: %s <input.fasta> <out_prefix> [k]\n"
                     "  reads <out_prefix>.unitigs.tsv\n"
                     "  writes <out_prefix>.occ.tsv and <out_prefix>.refs.tsv\n"
                     "  k defaults to 31; max %d.\n",
                     argv[0], kMaxK);
        return 1;
    }
    const std::string input_path = argv[1];
    const std::string out_prefix = argv[2];
    const int k = (argc == 4) ? std::atoi(argv[3]) : 31;
    if (k < 2 || k > kMaxK) {
        std::fprintf(stderr, "[build_occ] error: k must be in [2, %d]\n", kMaxK);
        return 1;
    }
    const std::uint64_t mask = ((static_cast<std::uint64_t>(1) << (2 * k)) - 1);

    const std::string unitigs_path = out_prefix + ".unitigs.tsv";
    ParsedUnitigs unitigs;
    try {
        unitigs = read_unitigs_tsv(unitigs_path);
    } catch (const std::exception& e) {
        std::fprintf(stderr, "[build_occ] %s\n", e.what());
        return 1;
    }
    std::fprintf(stderr, "[build_occ] read %zu unitig(s) from %s\n",
                 unitigs.seqs.size(), unitigs_path.c_str());

    std::unordered_map<std::uint64_t, KmerLoc> kmer_loc;
    kmer_loc.reserve(unitigs.seqs.size() * 4);
    std::size_t total_unitig_kmers = 0;
    for (std::size_t uid = 0; uid < unitigs.seqs.size(); ++uid) {
        const std::string& seq = unitigs.seqs[uid];
        if (static_cast<int>(seq.size()) < k) {
            continue;
        }
        std::uint64_t fwd = 0;
        std::uint64_t rev = 0;
        int filled = 0;
        for (int i = 0; i < static_cast<int>(seq.size()); ++i) {
            const int b = base_to_2bit(seq[i]);
            if (b < 0) {
                std::fprintf(stderr,
                             "[build_occ] error: unitig %zu has non-ACGT base at offset %d\n",
                             uid, i);
                return 1;
            }
            fwd = ((fwd << 2) | static_cast<std::uint64_t>(b)) & mask;
            rev = (rev >> 2) |
                  (static_cast<std::uint64_t>(b ^ 3) << (2 * (k - 1)));
            if (++filled >= k) {
                const std::uint64_t canon = std::min(fwd, rev);
                const char o_local = (fwd == canon) ? '+' : '-';
                const int off_local = i - (k - 1);
                const auto [it, inserted] = kmer_loc.try_emplace(
                    canon, KmerLoc{static_cast<int>(uid), off_local, o_local});
                if (!inserted) {
                    std::fprintf(stderr,
                                 "[build_occ] error: canonical k-mer appears in two unitigs"
                                 " (unitigs %d and %zu); the canonical compacted dBG invariant is broken\n",
                                 it->second.unitig_id, uid);
                    return 1;
                }
                ++total_unitig_kmers;
            }
        }
    }
    std::fprintf(stderr,
                 "[build_occ] indexed %zu canonical k-mers across all unitigs\n",
                 total_unitig_kmers);

    std::vector<std::string> ref_names;
    std::vector<std::string> ref_seqs;
    try {
        ref_seqs = read_fasta_with_names(input_path, ref_names);
    } catch (const std::exception& e) {
        std::fprintf(stderr, "[build_occ] %s\n", e.what());
        return 1;
    }
    std::fprintf(stderr,
                 "[build_occ] read %zu reference sequence(s) from %s\n",
                 ref_seqs.size(), input_path.c_str());

    std::vector<OccRow> rows;
    std::size_t missing_lookups = 0;

    for (std::size_t ref_id = 0; ref_id < ref_seqs.size(); ++ref_id) {
        const std::string& seq = ref_seqs[ref_id];
        std::uint64_t fwd = 0;
        std::uint64_t rev = 0;
        int filled = 0;

        bool have_run = false;
        OccRow run{};
        int run_last_off = 0;

        const auto flush_run = [&]() {
            if (have_run) {
                rows.push_back(run);
                have_run = false;
            }
        };

        for (int i = 0; i < static_cast<int>(seq.size()); ++i) {
            const int b = base_to_2bit(seq[i]);
            if (b < 0) {
                flush_run();
                filled = 0;
                fwd = 0;
                rev = 0;
                continue;
            }
            fwd = ((fwd << 2) | static_cast<std::uint64_t>(b)) & mask;
            rev = (rev >> 2) |
                  (static_cast<std::uint64_t>(b ^ 3) << (2 * (k - 1)));
            if (++filled < k) {
                continue;
            }
            const long long p = static_cast<long long>(i) - (k - 1);
            const std::uint64_t canon = std::min(fwd, rev);
            const auto it = kmer_loc.find(canon);
            if (it == kmer_loc.end()) {
                ++missing_lookups;
                flush_run();
                continue;
            }
            const int u_id = it->second.unitig_id;
            const int off_local = it->second.off_local;
            const char o_local = it->second.o_local;
            const char ref_strand = (fwd == canon) ? '+' : '-';
            const char o_ref = (o_local == ref_strand) ? '+' : '-';

            if (have_run && run.unitig_id == u_id && run.orient == o_ref) {
                const int expected = (o_ref == '+') ? run_last_off + 1
                                                    : run_last_off - 1;
                if (off_local == expected) {
                    run_last_off = off_local;
                    ++run.walk_len;
                    continue;
                }
            }

            flush_run();
            run = OccRow{u_id, static_cast<int>(ref_id), p, off_local, o_ref, 1};
            run_last_off = off_local;
            have_run = true;
        }
        flush_run();
    }

    if (missing_lookups > 0) {
        std::fprintf(stderr,
                     "[build_occ] WARNING: %zu reference k-mer(s) were not found in any unitig"
                     "; this means the FASTA passed to build_occ does not match the one used"
                     " to build the unitigs\n",
                     missing_lookups);
    }
    std::fprintf(stderr, "[build_occ] built %zu occurrence rows\n", rows.size());

    {
        long long ref_kmer_windows = 0;
        for (const auto& s : ref_seqs) {
            int filled = 0;
            for (char c : s) {
                if (base_to_2bit(c) < 0) {
                    filled = 0;
                } else if (++filled >= k) {
                    ++ref_kmer_windows;
                }
            }
        }
        long long covered = 0;
        for (const auto& r : rows) {
            covered += static_cast<long long>(r.walk_len);
        }
        if (covered != ref_kmer_windows) {
            std::fprintf(stderr,
                         "[build_occ] WARNING: covered=%lld vs ref_kmer_windows=%lld;"
                         " mismatch indicates a bug in the occurrence builder\n",
                         covered, ref_kmer_windows);
        } else {
            std::fprintf(stderr,
                         "[build_occ] sanity OK: %lld k-mer windows covered exactly by occurrence rows\n",
                         covered);
        }
    }

    {
        std::size_t mismatches = 0;
        for (const auto& r : rows) {
            const std::string& useq = unitigs.seqs[static_cast<std::size_t>(r.unitig_id)];
            const std::string& rseq = ref_seqs[static_cast<std::size_t>(r.ref_id)];
            const int L = static_cast<int>(useq.size());
            bool ok = true;
            for (int s = 0; s < r.walk_len; ++s) {
                const int u_off = (r.orient == '+') ? r.entry_off + s
                                                    : r.entry_off - s;
                const long long r_pos = r.ref_pos + s;
                if (u_off < 0 || u_off + k > L) { ok = false; break; }
                if (r_pos < 0 || r_pos + k > static_cast<long long>(rseq.size())) {
                    ok = false; break;
                }
                std::string u_kmer = useq.substr(static_cast<std::size_t>(u_off),
                                                 static_cast<std::size_t>(k));
                if (r.orient == '-') {
                    u_kmer = reverse_complement_string(u_kmer);
                }
                const std::string r_kmer = rseq.substr(static_cast<std::size_t>(r_pos),
                                                       static_cast<std::size_t>(k));
                if (u_kmer != r_kmer) { ok = false; break; }
            }
            if (!ok) {
                ++mismatches;
            }
        }
        if (mismatches > 0) {
            std::fprintf(stderr,
                         "[build_occ] WARNING: %zu / %zu rows do not reconstruct the FASTA exactly\n",
                         mismatches, rows.size());
        } else {
            std::fprintf(stderr,
                         "[build_occ] reconstruction OK: every k-mer in every run matches the FASTA\n");
        }
    }

    {
        const std::filesystem::path out_path(out_prefix);
        const auto parent = out_path.parent_path();
        if (!parent.empty()) {
            std::error_code ec;
            std::filesystem::create_directories(parent, ec);
            if (ec) {
                std::fprintf(stderr,
                             "[build_occ] error: cannot create output directory %s: %s\n",
                             parent.string().c_str(), ec.message().c_str());
                return 1;
            }
        }
    }

    const std::string occ_path = out_prefix + ".occ.tsv";
    const std::string refs_path = out_prefix + ".refs.tsv";
    std::ofstream o_out(occ_path);
    std::ofstream r_out(refs_path);
    if (!o_out.is_open() || !r_out.is_open()) {
        std::fprintf(stderr, "[build_occ] error: cannot write outputs under %s\n",
                     out_prefix.c_str());
        return 1;
    }

    o_out << "unitig_id\tref_id\tref_pos\tentry_off\torient\twalk_len\n";
    for (const auto& r : rows) {
        o_out << r.unitig_id << '\t'
              << r.ref_id << '\t'
              << r.ref_pos << '\t'
              << r.entry_off << '\t'
              << r.orient << '\t'
              << r.walk_len << '\n';
    }

    r_out << "ref_id\tname\tlength\n";
    for (std::size_t i = 0; i < ref_seqs.size(); ++i) {
        r_out << i << '\t'
              << (ref_names[i].empty() ? "<unnamed>" : ref_names[i]) << '\t'
              << ref_seqs[i].size() << '\n';
    }

    std::fprintf(stderr,
                 "[build_occ] wrote %s\n"
                 "[build_occ] wrote %s\n",
                 occ_path.c_str(), refs_path.c_str());
    return 0;
}
