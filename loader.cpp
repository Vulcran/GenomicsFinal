#include "loader.hpp"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <unordered_map>

namespace flat_index {

// Parse CIGAR string like "30M" → 30.
static uint8_t parse_cigar_overlap(const std::string& cigar) {
    if (cigar.empty() || cigar == "*") return 0;
    size_t end = 0;
    long v = std::stol(cigar, &end);
    if (end == 0 || end >= cigar.size() || cigar[end] != 'M')
        throw std::runtime_error("unrecognized CIGAR token: " + cigar);
    return static_cast<uint8_t>(v);
}

// Strip trailing '\r' for Windows-style line endings.
static void strip_cr(std::string& s) {
    if (!s.empty() && s.back() == '\r') s.pop_back();
}

// Parse occ.tsv into InOccurrence rows. Column order from build_occ:
//   unitig_id  ref_id  ref_pos  entry_off  orient  walk_len
static void parse_occ_stream(std::istream& in,
                              std::vector<InOccurrence>& occs) {
    std::string line;
    // Parse header to map column names → indices.
    if (!std::getline(in, line)) return;
    strip_cr(line);

    std::unordered_map<std::string, int> col;
    {
        std::istringstream hs(line);
        std::string tok;
        int idx = 0;
        while (std::getline(hs, tok, '\t')) col[tok] = idx++;
    }
    // Required column indices (with defaults for the known format).
    auto ci = [&](const std::string& name, int def) {
        auto it = col.find(name);
        return it != col.end() ? it->second : def;
    };
    int c_uid  = ci("unitig_id", 0);
    int c_rid  = ci("ref_id",    1);
    int c_rpos = ci("ref_pos",   2);
    int c_eoff = ci("entry_off", 3);
    int c_ori  = ci("orient",    4);
    int c_wlen = ci("walk_len",  5);
    int max_col = std::max({c_uid, c_rid, c_rpos, c_eoff, c_ori, c_wlen});

    while (std::getline(in, line)) {
        strip_cr(line);
        if (line.empty()) continue;
        std::vector<std::string> f;
        {
            std::istringstream ss(line);
            std::string tok;
            while (std::getline(ss, tok, '\t')) f.push_back(tok);
        }
        if (static_cast<int>(f.size()) <= max_col) continue;
        InOccurrence o;
        o.unitig_id = static_cast<uint32_t>(std::stoul(f[c_uid]));
        o.ref_id    = static_cast<uint16_t>(std::stoul(f[c_rid]));
        o.ref_pos   = static_cast<uint32_t>(std::stoul(f[c_rpos]));
        o.entry_off = static_cast<uint32_t>(std::stoul(f[c_eoff]));
        o.orient    = f[c_ori].empty() ? '+' : f[c_ori][0];
        o.walk_len  = static_cast<uint32_t>(std::stoul(f[c_wlen]));
        occs.push_back(o);
    }
}

LoadedData load_from_tsv(const std::string& pfx, bool require_occ) {
    LoadedData d;

    // --- unitigs.tsv: id  length  sequence ---
    {
        std::string path = pfx + ".unitigs.tsv";
        std::ifstream f(path);
        if (!f) throw std::runtime_error("cannot open: " + path);
        std::string line;
        std::getline(f, line); // skip header
        while (std::getline(f, line)) {
            strip_cr(line);
            if (line.empty()) continue;
            std::istringstream ss(line);
            std::string id_s, len_s, seq_s;
            if (!std::getline(ss, id_s, '\t')) continue;
            if (!std::getline(ss, len_s, '\t')) continue;
            if (!std::getline(ss, seq_s, '\t')) continue;
            d.unitigs.push_back({static_cast<uint32_t>(std::stoul(id_s)), std::move(seq_s)});
        }
        std::sort(d.unitigs.begin(), d.unitigs.end(),
                  [](const InUnitig& a, const InUnitig& b) { return a.id < b.id; });
    }

    // --- edges.tsv: from_id  from_orient  to_id  to_orient  overlap ---
    {
        std::string path = pfx + ".edges.tsv";
        std::ifstream f(path);
        if (!f) throw std::runtime_error("cannot open: " + path);
        std::string line;
        std::getline(f, line); // skip header
        while (std::getline(f, line)) {
            strip_cr(line);
            if (line.empty()) continue;
            std::istringstream ss(line);
            std::string fid, fo, tid, to_, ov;
            if (!std::getline(ss, fid, '\t')) continue;
            if (!std::getline(ss, fo,  '\t')) continue;
            if (!std::getline(ss, tid, '\t')) continue;
            if (!std::getline(ss, to_, '\t')) continue;
            if (!std::getline(ss, ov,  '\t')) continue;
            d.edges.push_back({
                static_cast<uint32_t>(std::stoul(fid)), fo.empty() ? '+' : fo[0],
                static_cast<uint32_t>(std::stoul(tid)), to_.empty() ? '+' : to_[0],
                parse_cigar_overlap(ov)
            });
        }
    }

    // --- occ.tsv (optional) ---
    {
        std::string path = pfx + ".occ.tsv";
        std::ifstream f(path);
        if (!f && require_occ) throw std::runtime_error("cannot open: " + path);
        if (f) parse_occ_stream(f, d.occs);
    }

    return d;
}

void load_occ_tsv(LoadedData& d, const std::string& path) {
    std::ifstream f(path);
    if (!f) throw std::runtime_error("cannot open: " + path);
    parse_occ_stream(f, d.occs);
}

LoadedData load_from_gfa(const std::string& gfa_path) {
    LoadedData d;
    std::ifstream f(gfa_path);
    if (!f) throw std::runtime_error("cannot open GFA: " + gfa_path);

    // Map GFA segment name → sequential uint32_t id.
    std::unordered_map<std::string, uint32_t> seg_id;
    uint32_t next_id = 0;

    auto get_id = [&](const std::string& name) -> uint32_t {
        auto [it, inserted] = seg_id.emplace(name, next_id);
        if (inserted) ++next_id;
        return it->second;
    };

    // Buffer edges until all segments are seen.
    struct RawEdge { std::string from, to; char fo, to_o; uint8_t ov; };
    std::vector<RawEdge> raw_edges;

    std::string line;
    while (std::getline(f, line)) {
        strip_cr(line);
        if (line.empty() || line[0] == '#') continue;

        if (line[0] == 'S') {
            // S  <name>  <sequence>  [tags...]
            std::istringstream ss(line);
            std::string tag, name, seq;
            ss >> tag >> name >> seq;
            if (seq == "*") continue;
            uint32_t id = get_id(name);
            // Grow unitig vector to hold this id.
            if (id >= d.unitigs.size()) d.unitigs.resize(id + 1);
            d.unitigs[id] = {id, std::move(seq)};

        } else if (line[0] == 'L') {
            // L  <from>  <from_orient>  <to>  <to_orient>  <CIGAR>
            std::istringstream ss(line);
            std::string tag, from, from_o, to, to_o, cigar;
            ss >> tag >> from >> from_o >> to >> to_o >> cigar;
            raw_edges.push_back({from, to,
                from_o.empty() ? '+' : from_o[0],
                to_o.empty()   ? '+' : to_o[0],
                parse_cigar_overlap(cigar)});
        }
        // H and P lines are ignored.
    }

    // Now that all segment IDs are assigned, build edges.
    d.edges.reserve(raw_edges.size());
    for (const auto& e : raw_edges) {
        uint32_t fid = get_id(e.from);
        uint32_t tid = get_id(e.to);
        d.edges.push_back({fid, e.fo, tid, e.to_o, e.ov});
    }

    return d;
}

} // namespace flat_index
