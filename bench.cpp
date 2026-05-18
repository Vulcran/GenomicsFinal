#include <iostream>
#include <vector>
#include <string>
#include <chrono>
#include <random>
#include <unordered_map>
#include "data_model.hpp"
#include "query.hpp"
#include "kmer.hpp"
#include "metrics.hpp"
#include "loader.hpp"

using namespace flat_index;

struct FlatListOccurrence {
    uint16_t ref_id;
    uint32_t ref_pos;
    bool orient;
};

template<typename K, typename V>
struct InstrumentedMap {
    struct Node {
        K key;
        V value;
        Node* next;
    };
    std::vector<Node*> buckets;
    size_t size_ = 0;
    
    InstrumentedMap(size_t n) : buckets(n, nullptr) {}
    ~InstrumentedMap() {
        for (auto b : buckets) {
            Node* curr = b;
            while (curr) {
                Node* next = curr->next;
                delete curr;
                curr = next;
            }
        }
    }
    
    void reserve(size_t n) {
        if (n <= buckets.size()) return;
        // Simple resize: just clear and allocate more buckets.
        // For benchmark, we know size in advance.
        buckets.assign(n, nullptr);
    }
    
    void add(K k, const typename V::value_type& v) {
        size_t h = std::hash<K>{}(k) % buckets.size();
        Node* curr = buckets[h];
        while (curr) {
            if (curr->key == k) {
                curr->value.push_back(v);
                return;
            }
            curr = curr->next;
        }
        buckets[h] = new Node{k, {v}, buckets[h]};
        size_++;
    }
    
    const V* find(K k) const {
        size_t h = std::hash<K>{}(k) % buckets.size();
        g_cache_counter.record_miss(); // Bucket access
        Node* curr = buckets[h];
        while (curr) {
            g_cache_counter.record_miss(); // Node access
            if (curr->key == k) {
                return &curr->value;
            }
            curr = curr->next;
        }
        return nullptr;
    }
    
    size_t num_buckets() const { return buckets.size(); }
    size_t size() const { return size_; }
};

std::string generate_random_dna(size_t len, std::mt19937& rng) {
    const char bases[] = "ACGT";
    std::string res;
    res.reserve(len);
    std::uniform_int_distribution<int> dist(0, 3);
    for (size_t i = 0; i < len; ++i) {
        res += bases[dist(rng)];
    }
    return res;
}

void run_benchmark(const std::string& tsv_prefix = "", int tsv_k = 31) {
    int k = tsv_k;
    std::vector<InUnitig> unitigs;
    std::vector<InOccurrence> occs;
    std::vector<InEdge> edges;

    if (!tsv_prefix.empty()) {
        std::cout << "Loading real data from: " << tsv_prefix << std::endl;
        auto data = load_from_tsv(tsv_prefix);
        unitigs = std::move(data.unitigs);
        occs    = std::move(data.occs);
        edges   = std::move(data.edges);
        std::cout << "Loaded " << unitigs.size() << " unitigs, "
                  << occs.size() << " occurrences, "
                  << edges.size() << " edges." << std::endl;
    } else {
        size_t num_unitigs = 10000;
        size_t unitig_len = 10000;
        std::mt19937 gen(42);
        std::cout << "Generating " << num_unitigs << " random unitigs..." << std::endl;
        for (size_t i = 0; i < num_unitigs; ++i) {
            if (i % 100 == 0 || i == num_unitigs - 1) {
                std::cout << "\rGenerating: " << (i + 1) << " / " << num_unitigs << std::flush;
            }
            std::string seq = generate_random_dna(unitig_len, gen);
            unitigs.push_back({(uint32_t)i, seq});
            occs.push_back({(uint32_t)i, (uint16_t)i, 0, 0, (uint8_t)(unitig_len - k + 1), '+'});
        }
        std::cout << std::endl;
    }
    
    std::cout << "Building Pufferfish Index..." << std::endl;
    auto t0 = std::chrono::steady_clock::now();
    FullIndex idx = build_full_index(unitigs, occs, edges, k);
    auto t1 = std::chrono::steady_clock::now();
    double p_build_time = std::chrono::duration<double>(t1 - t0).count();
    
    std::cout << "Building Flat List Index..." << std::endl;
    t0 = std::chrono::steady_clock::now();
    size_t total_kmers = 0;
    for (const auto& u : unitigs) total_kmers += (u.seq.length() - k + 1);
    
    InstrumentedMap<uint64_t, std::vector<FlatListOccurrence>> flat_list_map(total_kmers);
    for (size_t u_idx = 0; u_idx < unitigs.size(); ++u_idx) {
        if (u_idx % 100 == 0 || u_idx == unitigs.size() - 1) {
            std::cout << "\rIndexing: " << (u_idx + 1) << " / " << unitigs.size() << std::flush;
        }
        const auto& u = unitigs[u_idx];
        for (size_t i = 0; i + k <= u.seq.length(); ++i) {
            uint64_t kc = canonical(encode(u.seq.substr(i, k)), k);
            flat_list_map.add(kc, {(uint16_t)u.id, (uint32_t)i, false});
        }
    }
    std::cout << std::endl;
    t1 = std::chrono::steady_clock::now();
    double n_build_time = std::chrono::duration<double>(t1 - t0).count();

    // Collect all queryable (unitig_idx, max_pos) pairs so real data with
    // variable unitig lengths works correctly.
    std::vector<std::pair<size_t, size_t>> queryable;
    for (size_t ui = 0; ui < unitigs.size(); ++ui) {
        if (unitigs[ui].seq.size() >= static_cast<size_t>(k))
            queryable.push_back({ui, unitigs[ui].seq.size() - k});
    }

    size_t num_queries = 100000;
    std::vector<std::string> queries;
    std::mt19937 rng(43);
    std::uniform_int_distribution<size_t> u_dist(0, queryable.size() - 1);

    for (size_t i = 0; i < num_queries; ++i) {
        auto [u_idx, max_pos] = queryable[u_dist(rng)];
        std::uniform_int_distribution<size_t> p_dist(0, max_pos);
        size_t pos = p_dist(rng);
        queries.push_back(unitigs[u_idx].seq.substr(pos, k));
    }
    
    std::cout << "Running Pufferfish queries..." << std::endl;
    g_cache_counter.reset();
    t0 = std::chrono::steady_clock::now();
    size_t p_hits = 0;
    for (const auto& q : queries) {
        auto hits = query(idx, q);
        p_hits += hits.size();
    }
    t1 = std::chrono::steady_clock::now();
    double p_query_time = std::chrono::duration<double>(t1 - t0).count();
    size_t p_misses = g_cache_counter.misses;

    std::cout << "Running Flat List queries..." << std::endl;
    g_cache_counter.reset();
    t0 = std::chrono::steady_clock::now();
    size_t n_hits = 0;
    for (const auto& q : queries) {
        uint64_t kc = canonical(encode(q), k);
        auto val = flat_list_map.find(kc);
        if (val) {
            g_cache_counter.record_miss(); // Vector data access
            n_hits += val->size();
        }
    }
    t1 = std::chrono::steady_clock::now();
    double n_query_time = std::chrono::duration<double>(t1 - t0).count();
    size_t n_misses = g_cache_counter.misses;

    size_t p_mem = idx.flat_idx.flat.size() + idx.pos_table.size() * 8 + idx.mphf.table_size() * 4;
    // Estimate for InstrumentedMap: buckets + nodes + vectors + hits
    size_t n_mem = flat_list_map.num_buckets() * 8 + flat_list_map.size() * (8 + 24 + 8) + n_hits * 8;

    std::cout << "\nBenchmark Results (Correctness Hits: " << n_hits << " vs " << p_hits << ")" << std::endl;
    const int w_metric = 25;
    const int w_val = 20;

    auto print_sep = [&]() {
        std::cout << "+" << std::string(w_metric + 2, '-') << "+" 
                  << std::string(w_val + 2, '-') << "+" 
                  << std::string(w_val + 2, '-') << "+" << std::endl;
    };

    print_sep();
    printf("| %-*s | %*s | %*s |\n", w_metric, "Metric", w_val, "Flat List", w_val, "Pufferfish");
    print_sep();
    printf("| %-*s | %*.2f | %*.2f |\n", w_metric, "Throughput (M ops/s)", w_val, (num_queries / n_query_time) / 1e6, w_val, (num_queries / p_query_time) / 1e6);
    printf("| %-*s | %*.2f | %*.2f |\n", w_metric, "Memory (MB)", w_val, n_mem / 1e6, w_val, p_mem / 1e6);
    printf("| %-*s | %*.2f | %*.2f |\n", w_metric, "Build Time (s)", w_val, n_build_time, w_val, p_build_time);
    printf("| %-*s | %*.2f | %*.2f |\n", w_metric, "Cache Misses / Query", w_val, (double)n_misses / num_queries, w_val, (double)p_misses / num_queries);
    print_sep();
}

int main(int argc, char** argv) {
    std::string tsv_prefix;
    int k = 31;
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--tsv" && i + 1 < argc) tsv_prefix = argv[++i];
        else if (arg == "--k" && i + 1 < argc) k = std::stoi(argv[++i]);
        else { std::cerr << "Usage: bench [--tsv <prefix>] [--k <k>]\n"; return 1; }
    }
    run_benchmark(tsv_prefix, k);
    return 0;
}
