#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <optional>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

class SparsePufferfishLikeIndex {
public:
    struct QueryResult {
        std::string kmer;
        std::uint16_t unitig_id = 0;
        std::uint64_t unitig_bit_offset = 0;
        std::uint8_t local_offset = 0;
        std::string matched_unitig;
    };

    void build(const std::vector<std::string>& unitigs, const std::size_t k, const std::size_t sample_rate = 4) {
        clear();
        k_ = k;
        sample_rate_ = std::max<std::size_t>(1, sample_rate);
        unitigs_ = unitigs;
        resetCacheStats();

        for (std::size_t unitig_id = 0; unitig_id < unitigs.size(); ++unitig_id) {
            unitig_map_.push_back(static_cast<std::uint64_t>(unitig_flat_.size()));
            unitig_flat_ += unitigs[unitig_id];
        }

        std::unordered_map<std::string, std::uint64_t> first_occurrence;
        for (std::size_t unitig_id = 0; unitig_id < unitigs.size(); ++unitig_id) {
            const auto& unitig = unitigs[unitig_id];
            if (unitig.size() < k) {
                continue;
            }

            for (std::size_t i = 0; i + k <= unitig.size(); ++i) {
                const std::string kmer = unitig.substr(i, k);
                if (first_occurrence.find(kmer) == first_occurrence.end()) {
                    first_occurrence.emplace(kmer, pack56_8(unitig_id, static_cast<std::uint8_t>(i)));
                }
            }
        }

        std::vector<std::string> unique_kmers;
        unique_kmers.reserve(first_occurrence.size());
        for (std::unordered_map<std::string, std::uint64_t>::const_iterator it = first_occurrence.begin();
             it != first_occurrence.end();
             ++it) {
            unique_kmers.push_back(it->first);
        }

        mphf_.build(unique_kmers);
        slot_kmers_.assign(unique_kmers.size(), "");
        sampled_slot_.assign(unique_kmers.size(), false);
        sampled_pos_.assign(unique_kmers.size(), 0);
        anchor_slot_.assign(unique_kmers.size(), 0);
        delta_to_anchor_.assign(unique_kmers.size(), 0);

        for (std::size_t i = 0; i < unique_kmers.size(); ++i) {
            const std::string& kmer = unique_kmers[i];
            const std::size_t slot = mphf_.lookup(kmer);
            slot_kmers_[slot] = kmer;
        }

        for (std::size_t i = 0; i < unique_kmers.size(); ++i) {
            const std::string& kmer = unique_kmers[i];
            const std::size_t slot = mphf_.lookup(kmer);
            const std::uint64_t packed_occurrence = first_occurrence[kmer];
            const std::pair<std::uint64_t, std::uint8_t> unpacked = unpack56_8(packed_occurrence);
            const std::size_t unitig_id = static_cast<std::size_t>(unpacked.first);
            const std::size_t local_offset = static_cast<std::size_t>(unpacked.second);

            if (local_offset % sample_rate_ == 0) {
                sampled_slot_[slot] = true;
                sampled_pos_[slot] = packed_occurrence;
                anchor_slot_[slot] = slot;
                delta_to_anchor_[slot] = 0;
                continue;
            }

            const std::size_t anchor_local = (local_offset / sample_rate_) * sample_rate_;
            if (unitig_id >= unitigs_.size() || anchor_local + k_ > unitigs_[unitig_id].size()) {
                sampled_slot_[slot] = true;
                sampled_pos_[slot] = packed_occurrence;
                anchor_slot_[slot] = slot;
                delta_to_anchor_[slot] = 0;
                continue;
            }

            const std::string anchor_kmer = unitigs_[unitig_id].substr(anchor_local, k_);
            const std::size_t sampled_anchor_slot = mphf_.lookup(anchor_kmer);
            if (sampled_anchor_slot >= slot_kmers_.size() ||
                slot_kmers_[sampled_anchor_slot] != anchor_kmer ||
                anchor_local % sample_rate_ != 0) {
                sampled_slot_[slot] = true;
                sampled_pos_[slot] = packed_occurrence;
                anchor_slot_[slot] = slot;
                delta_to_anchor_[slot] = 0;
                continue;
            }

            sampled_slot_[slot] = false;
            anchor_slot_[slot] = static_cast<std::uint32_t>(sampled_anchor_slot);
            delta_to_anchor_[slot] = static_cast<std::int8_t>(static_cast<int>(local_offset) - static_cast<int>(anchor_local));

            if (!sampled_slot_[sampled_anchor_slot]) {
                sampled_slot_[sampled_anchor_slot] = true;
                sampled_pos_[sampled_anchor_slot] = pack56_8(unitig_id, static_cast<std::uint8_t>(anchor_local));
                anchor_slot_[sampled_anchor_slot] = static_cast<std::uint32_t>(sampled_anchor_slot);
                delta_to_anchor_[sampled_anchor_slot] = 0;
            }
        }
    }

    std::vector<QueryResult> queryAll(const std::string& kmer) const {
        if (k_ == 0 || kmer.size() != k_) {
            return {};
        }

        touchCache(kRegionMPHF, 0);
        if (!mphf_.contains(kmer)) {
            return {};
        }

        touchCache(kRegionMPHF, 1);
        const std::size_t slot = mphf_.lookup(kmer);
        touchCache(kRegionSlotKmers, slot);
        if (slot >= slot_kmers_.size() || slot_kmers_[slot] != kmer) {
            return {};
        }

        std::size_t resolved_slot = slot;
        std::size_t resolved_local_offset = 0;
        std::size_t resolved_unitig_id = 0;

        touchCache(kRegionSampledFlag, slot);
        if (sampled_slot_[slot]) {
            touchCache(kRegionSampledPos, slot);
            const std::pair<std::uint64_t, std::uint8_t> unpacked = unpack56_8(sampled_pos_[slot]);
            resolved_unitig_id = static_cast<std::size_t>(unpacked.first);
            resolved_local_offset = static_cast<std::size_t>(unpacked.second);
        } else {
            touchCache(kRegionAnchorSlot, slot);
            resolved_slot = static_cast<std::size_t>(anchor_slot_[slot]);
            touchCache(kRegionSampledFlag, resolved_slot);
            if (resolved_slot >= sampled_pos_.size() || !sampled_slot_[resolved_slot]) {
                return {};
            }
            touchCache(kRegionSampledPos, resolved_slot);
            const std::pair<std::uint64_t, std::uint8_t> anchor_unpacked = unpack56_8(sampled_pos_[resolved_slot]);
            resolved_unitig_id = static_cast<std::size_t>(anchor_unpacked.first);
            const int anchor_local = static_cast<int>(anchor_unpacked.second);
            touchCache(kRegionDelta, slot);
            const int delta = static_cast<int>(delta_to_anchor_[slot]);
            const int adjusted = anchor_local + delta;
            if (adjusted < 0) {
                return {};
            }
            resolved_local_offset = static_cast<std::size_t>(adjusted);
        }

        if (resolved_unitig_id >= unitigs_.size() || resolved_unitig_id >= unitig_map_.size()) {
            return {};
        }

        touchCache(kRegionUnitigMap, resolved_unitig_id);
        const std::string& matched_unitig = unitigs_[resolved_unitig_id];
        if (resolved_local_offset + k_ > matched_unitig.size()) {
            return {};
        }
        if (matched_unitig.substr(resolved_local_offset, k_) != kmer) {
            return {};
        }

        std::vector<QueryResult> results;

        for (std::size_t unitig_id = 0; unitig_id < unitigs_.size(); ++unitig_id) {
            touchCache(kRegionUnitigs, unitig_id);
            const std::string& unitig = unitigs_[unitig_id];
            if (unitig.size() < k_) {
                continue;
            }
            for (std::size_t i = 0; i + k_ <= unitig.size(); ++i) {
                if (unitig.compare(i, k_, kmer) != 0) {
                    continue;
                }
                QueryResult result;
                result.kmer = kmer;
                result.unitig_id = static_cast<std::uint16_t>(unitig_id);
                touchCache(kRegionUnitigMap, unitig_id);
                result.unitig_bit_offset = unitig_map_[unitig_id];
                result.local_offset = static_cast<std::uint8_t>(i);
                result.matched_unitig = unitig;
                results.push_back(std::move(result));
            }
        }

        return results;
    }

    [[nodiscard]] std::size_t unitigFlatSize() const { return unitig_flat_.size(); }
    [[nodiscard]] std::size_t kmerFlatSize() const { return sampled_pos_.size(); }
    [[nodiscard]] std::size_t uniqueKmerCount() const { return slot_kmers_.size(); }
    [[nodiscard]] std::size_t sampledKmerCount() const {
        return static_cast<std::size_t>(std::count(sampled_slot_.begin(), sampled_slot_.end(), true));
    }

    struct CacheStats {
        std::size_t hits = 0;
        std::size_t misses = 0;
    };

    [[nodiscard]] CacheStats cacheStats() const {
        return {cache_hits_, cache_misses_};
    }

    void resetCacheStats() {
        cache_hits_ = 0;
        cache_misses_ = 0;
        std::fill(cache_valid_.begin(), cache_valid_.end(), false);
    }

private:
    class StaticMPHF {
    public:
        void build(const std::vector<std::string>& keys) {
            size_ = keys.size();
            if (size_ == 0) {
                seeds_.clear();
                slots_.clear();
                return;
            }

            seeds_.assign(size_, 0);
            slots_.assign(size_, std::numeric_limits<std::size_t>::max());

            std::vector<std::vector<std::size_t>> buckets(size_);
            for (std::size_t i = 0; i < keys.size(); ++i) {
                buckets[baseHash(keys[i]) % size_].push_back(i);
            }

            std::vector<std::size_t> order(size_);
            for (std::size_t i = 0; i < size_; ++i) {
                order[i] = i;
            }

            std::sort(order.begin(), order.end(), [&](const std::size_t a, const std::size_t b) {
                return buckets[a].size() > buckets[b].size();
            });

            std::vector<bool> used(size_, false);

            for (const auto bucket_idx : order) {
                const auto& bucket = buckets[bucket_idx];
                if (bucket.empty()) {
                    continue;
                }

                if (bucket.size() == 1) {
                    std::size_t slot = 0;
                    while (slot < size_ && used[slot]) {
                        ++slot;
                    }
                    used[slot] = true;
                    slots_[slot] = bucket[0];
                    seeds_[bucket_idx] = -static_cast<std::int32_t>(slot + 1);
                    continue;
                }

                std::int32_t seed = 1;
                std::vector<std::size_t> chosen;
                chosen.reserve(bucket.size());

                while (true) {
                    chosen.clear();
                    bool collision = false;

                    for (const auto key_index : bucket) {
                        const std::size_t slot = seededHash(keys[key_index], seed) % size_;
                        if (used[slot] || std::find(chosen.begin(), chosen.end(), slot) != chosen.end()) {
                            collision = true;
                            break;
                        }
                        chosen.push_back(slot);
                    }

                    if (!collision) {
                        for (std::size_t i = 0; i < bucket.size(); ++i) {
                            const auto slot = chosen[i];
                            used[slot] = true;
                            slots_[slot] = bucket[i];
                        }
                        seeds_[bucket_idx] = seed;
                        break;
                    }
                    ++seed;
                }
            }
        }

        [[nodiscard]] bool contains(const std::string& key) const {
            return size_ != 0 && !seeds_.empty() && slotForKey(key) < size_;
        }

        [[nodiscard]] std::size_t lookup(const std::string& key) const { return slotForKey(key); }

    private:
        static std::uint64_t baseHash(const std::string& s) {
            return std::hash<std::string>{}(s);
        }

        static std::uint64_t seededHash(const std::string& s, const std::int32_t seed) {
            return baseHash(s) ^ (0x9e3779b97f4a7c15ULL * static_cast<std::uint64_t>(seed));
        }

        [[nodiscard]] std::size_t slotForKey(const std::string& key) const {
            if (size_ == 0) {
                return size_;
            }
            const std::size_t bucket = baseHash(key) % size_;
            const std::int32_t seed = seeds_[bucket];
            if (seed < 0) {
                return static_cast<std::size_t>(-seed - 1);
            }
            return seededHash(key, seed) % size_;
        }

        std::size_t size_ = 0;
        std::vector<std::int32_t> seeds_;
        std::vector<std::size_t> slots_;
    };

    static std::uint64_t pack56_8(const std::uint64_t upper_56, const std::uint8_t lower_8) {
        const std::uint64_t mask_56 = (1ULL << 56U) - 1ULL;
        return ((upper_56 & mask_56) << 8U) | lower_8;
    }

    static std::pair<std::uint64_t, std::uint8_t> unpack56_8(const std::uint64_t packed) {
        return {packed >> 8U, static_cast<std::uint8_t>(packed & 0xFFU)};
    }

    void touchCache(const std::uint64_t region, const std::uint64_t index) const {
        const std::uint64_t address = (region << 48U) ^ index;
        const std::size_t line = static_cast<std::size_t>(address % kCacheLines);
        const std::uint64_t tag = address / kCacheLines;
        if (cache_valid_[line] && cache_tags_[line] == tag) {
            ++cache_hits_;
            return;
        }
        ++cache_misses_;
        cache_valid_[line] = true;
        cache_tags_[line] = tag;
    }

    void clear() {
        unitig_flat_.clear();
        unitig_map_.clear();
        sampled_slot_.clear();
        sampled_pos_.clear();
        anchor_slot_.clear();
        delta_to_anchor_.clear();
        slot_kmers_.clear();
        unitigs_.clear();
        k_ = 0;
        sample_rate_ = 1;
        mphf_ = StaticMPHF{};
    }

    std::size_t k_ = 0;
    std::size_t sample_rate_ = 1;
    std::vector<std::string> unitigs_;

    static constexpr std::size_t kCacheLines = 4096;
    static constexpr std::uint64_t kRegionMPHF = 1;
    static constexpr std::uint64_t kRegionSlotKmers = 2;
    static constexpr std::uint64_t kRegionSampledFlag = 3;
    static constexpr std::uint64_t kRegionSampledPos = 4;
    static constexpr std::uint64_t kRegionAnchorSlot = 5;
    static constexpr std::uint64_t kRegionDelta = 6;
    static constexpr std::uint64_t kRegionUnitigMap = 7;
    static constexpr std::uint64_t kRegionUnitigs = 8;

    mutable std::size_t cache_hits_ = 0;
    mutable std::size_t cache_misses_ = 0;
    mutable std::array<bool, kCacheLines> cache_valid_{};
    mutable std::array<std::uint64_t, kCacheLines> cache_tags_{};

    std::string unitig_flat_;
    std::vector<std::uint64_t> unitig_map_;

    std::vector<bool> sampled_slot_;
    std::vector<std::uint64_t> sampled_pos_;
    std::vector<std::uint32_t> anchor_slot_;
    std::vector<std::int8_t> delta_to_anchor_;
    std::vector<std::string> slot_kmers_;
    StaticMPHF mphf_;
};

namespace {
    struct GraphInput {
        std::size_t k = 0;
        std::vector<std::string> unitigs;
    };


    std::optional<GraphInput> loadCompactedColorGraph(const std::string& graph_path, std::ostream& err) {
        std::ifstream input(graph_path);
        if (!input) {
            err << "Unable to open graph file: " << graph_path << '\n';
            return std::nullopt;
        }

        GraphInput graph;
        std::string line;
        bool saw_k = false;
        while (std::getline(input, line)) {
            if (line.empty()) {
                continue;
            }

            if (!saw_k) {
                std::istringstream iss(line);
                if (!(iss >> graph.k) || graph.k == 0) {
                    err << "Invalid graph format: first non-empty line must be a positive k value\n";
                    return std::nullopt;
                }
                saw_k = true;
                continue;
            }

            graph.unitigs.push_back(line);
        }

        if (!saw_k) {
            err << "Invalid graph format: missing k value\n";
            return std::nullopt;
        }
        if (graph.unitigs.empty()) {
            err << "Invalid graph format: missing unitig lines\n";
            return std::nullopt;
        }

        return graph;
    }

    std::optional<std::vector<std::string>> loadQueries(const std::string& queries_path, std::ostream& err) {
        std::ifstream input(queries_path);
        if (!input) {
            err << "Unable to open query file: " << queries_path << '\n';
            return std::nullopt;
        }

        std::vector<std::string> queries;
        std::string line;
        while (std::getline(input, line)) {
            if (!line.empty()) {
                queries.push_back(line);
            }
        }

        if (queries.empty()) {
            err << "Query file contains no queries: " << queries_path << '\n';
            return std::nullopt;
        }

        return queries;
    }

    void writeQueryResults(std::ostream& out,
                           const std::string& query,
                           const std::vector<SparsePufferfishLikeIndex::QueryResult>& results) {
        if (results.empty()) {
            out << "Query " << query << " => NOT FOUND\n";
            return;
        }

        out << "Query " << query << " => FOUND " << results.size() << " hit(s)\n";
        for (std::size_t i = 0; i < results.size(); ++i) {
            const auto& result = results[i];
            out << "   hit " << (i + 1) << ":\n";
            out << "      unitig_id: " << result.unitig_id << '\n';
            out << "      unitig_char_offset: " << result.unitig_bit_offset << '\n';
            out << "      local_offset: " << static_cast<int>(result.local_offset) << '\n';
            out << "      unitig_prefix: "
                << result.matched_unitig.substr(0, std::min<std::size_t>(result.matched_unitig.size(), 40))
                << '\n';
        }
    }
}

int main(int argc, char* argv[]) {
    if (argc != 4) {
        std::cerr << "Usage: " << argv[0] << " <graph_input_file> <queries_input_file> <output_file>\n";
        return 1;
    }

    const std::string graph_input_path = argv[1];
    const std::string queries_input_path = argv[2];
    const std::string output_path = argv[3];

    const auto graph_input = loadCompactedColorGraph(graph_input_path, std::cerr);
    if (!graph_input.has_value()) {
        return 1;
    }

    const auto queries = loadQueries(queries_input_path, std::cerr);
    if (!queries.has_value()) {
        return 1;
    }

    SparsePufferfishLikeIndex index;
    index.build(graph_input->unitigs, graph_input->k, 4);

    std::ofstream output(output_path);
    if (!output) {
        std::cerr << "Unable to open output file: " << output_path << '\n';
        return 1;
    }

    output << "Sparse Pufferfish-like index built\n";
    output << " - unitigFlat chars: " << index.unitigFlatSize() << '\n';
    output << " - unique kmers: " << index.uniqueKmerCount() << '\n';
    output << " - sampled kmers: " << index.sampledKmerCount() << '\n';
    output << " - kmerFlat entries: " << index.kmerFlatSize() << "\n\n";

    for (const auto& query : *queries) {
        const auto results = index.queryAll(query);
        writeQueryResults(output, query, results);
    }

    const SparsePufferfishLikeIndex::CacheStats stats = index.cacheStats();

    output << "\nCache summary:\n";
    output << " - cache hits: " << stats.hits << '\n';
    output << " - cache misses: " << stats.misses << '\n';

    return 0;
}
