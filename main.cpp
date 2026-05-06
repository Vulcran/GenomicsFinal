#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <optional>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

class FlatArrayIndex {
public:
    struct QueryResult {
        std::string kmer;
        std::uint16_t unitig_id = 0;
        std::uint64_t unitig_bit_offset = 0;
        std::uint8_t local_offset = 0;
        std::string matched_unitig;
    };

    void build(const std::vector<std::string>& unitigs, const std::size_t k) {
        clear();
        k_ = k;
        unitigs_ = unitigs;

        for (std::size_t unitig_id = 0; unitig_id < unitigs.size(); ++unitig_id) {
            unitig_map_.push_back(static_cast<std::uint64_t>(unitig_flat_.size()));
            unitig_flat_ += unitigs[unitig_id];
        }

        std::unordered_map<std::string, std::vector<std::uint64_t>> occurrence_lists;
        for (std::size_t unitig_id = 0; unitig_id < unitigs.size(); ++unitig_id) {
            const auto& unitig = unitigs[unitig_id];
            if (unitig.size() < k) {
                continue;
            }

            for (std::size_t i = 0; i + k <= unitig.size(); ++i) {
                const std::string kmer = unitig.substr(i, k);
                const auto packed_occurrence = pack56_8(unitig_id, static_cast<std::uint8_t>(i));
                occurrence_lists[kmer].push_back(packed_occurrence);
            }
        }

        std::vector<std::string> unique_kmers;
        unique_kmers.reserve(occurrence_lists.size());
        for (const auto& [kmer, _] : occurrence_lists) {
            unique_kmers.push_back(kmer);
        }
        mphf_.build(unique_kmers);

        kmer_map_.assign(unique_kmers.size(), 0);
        slot_kmers_.assign(unique_kmers.size(), "");

        for (const auto& kmer : unique_kmers) {
            const std::size_t slot = mphf_.lookup(kmer);
            const auto& occurrences = occurrence_lists[kmer];
            const std::uint64_t start_index = static_cast<std::uint64_t>(kmer_flat_.size());
            const std::uint8_t count = static_cast<std::uint8_t>(std::min<std::size_t>(occurrences.size(), 255));
            kmer_map_[slot] = pack56_8(start_index, count);
            slot_kmers_[slot] = kmer;

            for (const auto packed_occurrence : occurrences) {
                kmer_flat_.push_back(packed_occurrence);
            }
        }
    }

    std::vector<QueryResult> queryAll(const std::string& kmer) const {
        if (!mphf_.contains(kmer)) {
            return {};
        }

        const std::size_t slot = mphf_.lookup(kmer);
        if (slot >= slot_kmers_.size() || slot_kmers_[slot] != kmer) {
            return {};
        }

        const auto [start_index, count] = unpack56_8(kmer_map_[slot]);
        const std::size_t start = static_cast<std::size_t>(start_index);
        const std::size_t end_index = start + count;
        if (end_index > kmer_flat_.size()) {
            return {};
        }

        std::vector<QueryResult> results;
        results.reserve(count);

        for (std::size_t i = start; i < end_index; ++i) {
            const auto [unitig_map_index, local_offset] = unpack56_8(kmer_flat_[i]);
            if (unitig_map_index >= unitig_map_.size()) {
                continue;
            }

            const std::uint16_t unitig_id = static_cast<std::uint16_t>(unitig_map_index);
            if (unitig_id >= unitigs_.size()) {
                continue;
            }

            QueryResult result;
            result.kmer = kmer;
            result.unitig_id = unitig_id;
            result.unitig_bit_offset = unitig_map_[unitig_id];
            result.local_offset = local_offset;
            result.matched_unitig = unitigs_[unitig_id];
            results.push_back(std::move(result));
        }

        return results;
    }

    void printUnitigStorage() const {
        std::cout << "unitigFlat: '" << unitig_flat_ << "'\n";
        std::cout << "unitigMap (unitig_id -> char offset):\n";
        for (std::size_t i = 0; i < unitig_map_.size(); ++i) {
            std::cout << "  [" << i << "] = " << unitig_map_[i] << "\n";
        }
        std::cout << '\n';
    }

    void printKmerStorage() const {
        std::cout << "kmerFlat (packed unitigMap_index/local_offset entries):\n";
        for (std::size_t i = 0; i < kmer_flat_.size(); ++i) {
            const auto [unitig_map_index, local_offset] = unpack56_8(kmer_flat_[i]);
            std::cout << "  [" << i << "] = (" << unitig_map_index << "," << static_cast<int>(local_offset) << ")\n";
        }
        std::cout << '\n';

        std::cout << "kmerMap (MPHF slot -> packed start_index/count):\n";
        for (std::size_t i = 0; i < kmer_map_.size(); ++i) {
            const auto [start_index, count] = unpack56_8(kmer_map_[i]);
            std::cout << "  slot " << i
                      << " | kmer='" << slot_kmers_[i] << "'"
                      << " | start_index=" << start_index
                      << " | count=" << static_cast<int>(count) << "\n";
        }
        std::cout << '\n';
    }

    std::size_t unitigFlatSize() const { return unitig_flat_.size(); }
    std::size_t kmerFlatSize() const { return kmer_flat_.size(); }
    std::size_t uniqueKmerCount() const { return kmer_map_.size(); }

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

        [[nodiscard]] std::size_t lookup(const std::string& key) const {
            return slotForKey(key);
        }

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

    void clear() {
        unitig_flat_.clear();
        unitig_map_.clear();
        kmer_flat_.clear();
        kmer_map_.clear();
        slot_kmers_.clear();
        unitigs_.clear();
        k_ = 0;
        mphf_ = StaticMPHF{};
    }

    std::size_t k_ = 0;
    std::vector<std::string> unitigs_;

    std::string unitig_flat_;
    std::vector<std::uint64_t> unitig_map_;

    std::vector<std::uint64_t> kmer_flat_;
    std::vector<std::uint64_t> kmer_map_;
    std::vector<std::string> slot_kmers_;
    StaticMPHF mphf_;
};

namespace {
    struct GraphInput {
        std::size_t k = 0;
        std::vector<std::string> unitigs;
    };

    struct CacheStats {
        std::size_t hits = 0;
        std::size_t misses = 0;
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
                           const std::vector<FlatArrayIndex::QueryResult>& results) {
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

    FlatArrayIndex index;
    index.build(graph_input->unitigs, graph_input->k);

    std::ofstream output(output_path);
    if (!output) {
        std::cerr << "Unable to open output file: " << output_path << '\n';
        return 1;
    }

    output << "Flat-array + MPHF index built\n";
    output << " - unitigFlat chars: " << index.unitigFlatSize() << '\n';
    output << " - unique kmers: " << index.uniqueKmerCount() << '\n';
    output << " - kmerFlat entries: " << index.kmerFlatSize() << "\n\n";

    CacheStats stats;
    std::unordered_map<std::string, std::vector<FlatArrayIndex::QueryResult>> query_cache;
    for (const auto& query : *queries) {
        const auto it = query_cache.find(query);
        if (it != query_cache.end()) {
            ++stats.hits;
            writeQueryResults(output, query, it->second);
            continue;
        }

        ++stats.misses;
        auto results = index.queryAll(query);
        writeQueryResults(output, query, results);
        query_cache.emplace(query, std::move(results));
    }

    output << "\nCache summary:\n";
    output << " - cache hits: " << stats.hits << '\n';
    output << " - cache misses: " << stats.misses << '\n';

    return 0;
}