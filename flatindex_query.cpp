// flatindex_query — load a serialized flat index and query k-mers or reads.
//
// Usage:
//   flatindex_query <index.flat> [options]
//
// Options:
//   --kmer  <KMER>      Query a single k-mer string
//   --read  <READ>      Align a read (seed on first k-mer, walk the graph)
//   --fastq <file.fq>   Query every read in a FASTQ file
//   --stdin             Read one k-mer or read per line from stdin
//
// 
// Example:
//   flatindex_query 200bp.flat --kmer ACGTACGTACGTACGTACGTACGTACGTACG

#include <iostream>
#include <fstream>
#include <string>
#include <stdexcept>
#include "query.hpp"
#include "serializer.hpp"
#include "metrics.hpp"

using namespace flat_index;

static void print_hits(const std::vector<RefHit>& hits, const std::string& kmer) {
    if (hits.empty()) {
        std::cout << "  " << kmer << ": no hits\n";
        return;
    }
    for (const auto& h : hits) {
        std::cout << "  " << kmer << "  ref=" << h.ref_id
                  << "  pos=" << h.ref_pos
                  << "  orient=" << (h.orient ? '-' : '+') << "\n";
    }
}

static void print_alns(const std::vector<Alignment>& alns, const std::string& read) {
    if (alns.empty()) {
        std::cout << "  " << read << ": no alignment\n";
        return;
    }
    for (const auto& a : alns) {
        std::cout << "  read_len=" << read.size()
                  << "  unitig=" << a.unitig_id
                  << "  start_kmer=" << a.start_kmer_idx
                  << "  len_kmers=" << a.len_kmers
                  << "  orient=" << (a.orient ? '-' : '+') << "\n";
    }
}

// Parse a single FASTQ record from the stream. Returns false at EOF.
static bool read_fastq_record(std::istream& in, std::string& name, std::string& seq) {
    std::string line;
    // '@' header line
    while (std::getline(in, line)) {
        if (!line.empty() && line[0] == '@') { name = line.substr(1); break; }
    }
    if (in.eof()) return false;
    // Sequence
    if (!std::getline(in, seq)) return false;
    // '+' separator + quality (skip)
    if (!std::getline(in, line)) return false;
    if (!std::getline(in, line)) return false;
    return !seq.empty();
}

static void usage(const char* prog) {
    std::cerr
        << "Usage: " << prog << " <index.flat> [--kmer K] [--read R] [--fastq F] [--stdin]\n";
}

int main(int argc, char** argv) {
    if (argc < 2) { usage(argv[0]); return 1; }
    std::string index_path = argv[1];

    std::string kmer_arg, read_arg, fastq_path, kmer_file;
    bool from_stdin = false;
    bool no_output = false;

    for (int i = 2; i < argc; ++i) {
        std::string arg = argv[i];
        auto next = [&]() -> std::string {
            if (i + 1 >= argc) throw std::runtime_error("missing argument after " + arg);
            return argv[++i];
        };
        if      (arg == "--kmer")       kmer_arg   = next();
        else if (arg == "--read")       read_arg   = next();
        else if (arg == "--fastq")      fastq_path = next();
        else if (arg == "--kmer-file")  kmer_file  = next();
        else if (arg == "--stdin")      from_stdin = true;
        else if (arg == "--noOutput")   no_output  = true;
        else { std::cerr << "unknown argument: " << arg << "\n"; usage(argv[0]); return 1; }
    }

    FullIndex idx;
    try {
        std::cerr << "Loading index: " << index_path << " ... " << std::flush;
        idx = IndexSerializer::load(index_path);
        std::cerr << "ok (" << idx.flat_idx.flat.size() << " bytes flat, k="
                  << idx.flat_idx.k << ")\n";
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << "\n";
        return 1;
    }

    try {
        if (!kmer_arg.empty()) {
            auto hits = query(idx, kmer_arg);
            print_hits(hits, kmer_arg);

        } else if (!read_arg.empty()) {
            auto alns = align_read(idx, read_arg);
            print_alns(alns, read_arg);

        } else if (!kmer_file.empty()) {
            std::ifstream kf(kmer_file);
            if (!kf) throw std::runtime_error("cannot open: " + kmer_file);
            std::string kmer;
            size_t n = 0, hits = 0;
            while (std::getline(kf, kmer)) {
                if (kmer.empty()) continue;
                auto h = query(idx, kmer);
                if (!no_output) print_hits(h, kmer);
                if (!h.empty()) ++hits;
                ++n;
            }
            std::cerr << "queried " << n << " k-mers, " << hits << " hits\n";

        } else if (!fastq_path.empty()) {
            std::ifstream fq(fastq_path);
            if (!fq) throw std::runtime_error("cannot open: " + fastq_path);
            std::string name, seq;
            size_t n = 0;
            while (read_fastq_record(fq, name, seq)) {
                auto alns = align_read(idx, seq);
                if (!no_output) {
                    std::cout << "@" << name << "\n";
                    print_alns(alns, seq);
                }
                ++n;
            }
            std::cerr << "aligned " << n << " reads\n";
            std::cerr << "cache misses (simulated): " << g_cache_counter.misses
                      << "  (" << (g_cache_counter.misses / (n ? n : 1)) << " per read)\n";

        } else {
            // Interactive / --stdin: one sequence per line.
            std::istream& in = from_stdin ? std::cin : std::cin;
            std::cout << "Enter k-mers or reads (one per line, Ctrl-D to quit):\n";
            std::string line;
            while (std::getline(in, line)) {
                if (line.empty()) continue;
                int k = idx.flat_idx.k;
                if (static_cast<int>(line.size()) == k) {
                    auto hits = query(idx, line);
                    print_hits(hits, line);
                } else if (static_cast<int>(line.size()) > k) {
                    auto alns = align_read(idx, line);
                    print_alns(alns, line);
                } else {
                    std::cout << "  (sequence shorter than k=" << k << ", skipped)\n";
                }
            }
        }
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << "\n";
        return 1;
    }
    return 0;
}
