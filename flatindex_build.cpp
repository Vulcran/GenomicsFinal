// flatindex_build — build a flat index from TSV or GFA input and serialize it.
//
// Usage:
//   # From Project/ TSV pipeline (build_unitigs + build_occ output):
//   flatindex_build --tsv <prefix> --k <k> --out <index.flat>
//
//   # From a Cuttlefish GFA (unitigs + edges only; add --occ for occurrences):
//   flatindex_build --gfa <file.gfa> [--occ <occ.tsv>] --k <k> --out <index.flat>
//
// Example:
//   flatindex_build --tsv ../Project/data/200bp --k 31 --out 200bp.flat

#include <iostream>
#include <string>
#include <stdexcept>
#include "loader.hpp"
#include "query.hpp"
#include "serializer.hpp"

using namespace flat_index;

static void usage(const char* prog) {
    std::cerr
        << "Usage:\n"
        << "  " << prog << " --tsv <prefix> --k <k> --out <output.flat>\n"
        << "  " << prog << " --gfa <file.gfa> [--occ <occ.tsv>] --k <k> --out <output.flat>\n";
}

int main(int argc, char** argv) {
    std::string tsv_prefix, gfa_path, occ_path, out_path;
    int k = 0;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        auto next = [&]() -> std::string {
            if (i + 1 >= argc) throw std::runtime_error("missing argument after " + arg);
            return argv[++i];
        };
        if      (arg == "--tsv") tsv_prefix = next();
        else if (arg == "--gfa") gfa_path   = next();
        else if (arg == "--occ") occ_path   = next();
        else if (arg == "--k")   k          = std::stoi(next());
        else if (arg == "--out") out_path   = next();
        else { std::cerr << "unknown argument: " << arg << "\n"; usage(argv[0]); return 1; }
    }

    if (out_path.empty() || k <= 0) {
        std::cerr << "error: --out and --k are required\n";
        usage(argv[0]);
        return 1;
    }
    if (tsv_prefix.empty() == gfa_path.empty()) {
        std::cerr << "error: specify exactly one of --tsv or --gfa\n";
        usage(argv[0]);
        return 1;
    }

    try {
        LoadedData data;
        if (!tsv_prefix.empty()) {
            std::cout << "Loading TSV from prefix: " << tsv_prefix << std::endl;
            data = load_from_tsv(tsv_prefix);
        } else {
            std::cout << "Loading GFA: " << gfa_path << std::endl;
            data = load_from_gfa(gfa_path);
            if (!occ_path.empty()) {
                std::cout << "Loading occurrences: " << occ_path << std::endl;
                load_occ_tsv(data, occ_path);
            }
        }

        std::cout << "Loaded: " << data.unitigs.size() << " unitigs, "
                  << data.occs.size()   << " occurrences, "
                  << data.edges.size()  << " edges" << std::endl;

        std::cout << "Building index (k=" << k << ")..." << std::endl;
        FullIndex idx = build_full_index(data.unitigs, data.occs, data.edges, k);

        std::cout << "Saving to: " << out_path << std::endl;
        IndexSerializer::save(idx, out_path);

        std::cout << "Done. Flat array: "
                  << idx.flat_idx.flat.size() << " bytes, "
                  << idx.pos_table.size() << " POS slots." << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << "\n";
        return 1;
    }
    return 0;
}
