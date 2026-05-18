#ifndef LOADER_HPP
#define LOADER_HPP

#include "data_model.hpp"
#include <string>
#include <vector>

namespace flat_index {

struct LoadedData {
    std::vector<InUnitig>     unitigs;
    std::vector<InOccurrence> occs;
    std::vector<InEdge>       edges;
};

// Load from the four TSV files produced by build_unitigs + build_occ.
// path_prefix: e.g. "data/200bp" reads "data/200bp.unitigs.tsv", etc.
// If require_occ is false, missing occ.tsv is silently ignored.
LoadedData load_from_tsv(const std::string& path_prefix, bool require_occ = true);

// Load unitigs + edges from a GFA1 file (Cuttlefish --output-format gfa or similar).
// Occurrences are not populated; run build_occ separately and merge via load_occ_tsv.
LoadedData load_from_gfa(const std::string& gfa_path);

// Load only occurrences from a standalone occ.tsv into an existing LoadedData.
void load_occ_tsv(LoadedData& d, const std::string& occ_path);

} // namespace flat_index

#endif // LOADER_HPP
