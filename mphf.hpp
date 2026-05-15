#ifndef MPHF_HPP
#define MPHF_HPP

#include <vector>
#include <string>
#include <algorithm>
#include <limits>
#include <cstdint>
#include <iostream>
#include "metrics.hpp"

namespace flat_index {

/**
 * Minimal Perfect Hash Function.
 * Adapted from the original main.cpp to work with uint64_t keys (canonical k-mers).
 */
class StaticMPHF {
public:
    void build(const std::vector<uint64_t>& keys) {
        size_ = keys.size() * 1.2; // 20% slack
        if (keys.empty()) {
            size_ = 0;
            seeds_.clear();
            return;
        }

        seeds_.assign(size_, 0);
        
        std::cout << "\rMPHF: Bucketing keys... " << std::flush;
        std::vector<size_t> bucket_offsets(size_ + 1, 0);
        for (uint64_t k : keys) {
            bucket_offsets[baseHash(k) % size_ + 1]++;
        }
        for (size_t i = 0; i < size_; ++i) {
            bucket_offsets[i + 1] += bucket_offsets[i];
        }

        std::vector<size_t> bucket_data(keys.size());
        std::vector<size_t> current_pos = bucket_offsets;
        for (size_t i = 0; i < keys.size(); ++i) {
            size_t h = baseHash(keys[i]) % size_;
            bucket_data[current_pos[h]++] = i;
        }

        std::cout << "\rMPHF: Sorting buckets...  " << std::flush;
        std::vector<std::size_t> order(size_);
        for (std::size_t i = 0; i < size_; ++i) {
            order[i] = i;
        }

        std::sort(order.begin(), order.end(), [&](const std::size_t a, const std::size_t b) {
            size_t size_a = bucket_offsets[a + 1] - bucket_offsets[a];
            size_t size_b = bucket_offsets[b + 1] - bucket_offsets[b];
            return size_a > size_b;
        });

        std::vector<bool> used(size_, false);
        size_t next_free_slot = 0;

        size_t processed_buckets = 0;
        for (const auto bucket_idx : order) {
            if (processed_buckets % 10000 == 0 || processed_buckets == size_ - 1) {
                std::cout << "\rMPHF: Finding seeds... " << (processed_buckets + 1) << " / " << size_ << std::flush;
            }
            processed_buckets++;

            size_t b_start = bucket_offsets[bucket_idx];
            size_t b_end = bucket_offsets[bucket_idx + 1];
            size_t b_size = b_end - b_start;

            if (b_size == 0) continue;

            if (b_size == 1) {
                while (next_free_slot < size_ && used[next_free_slot]) {
                    ++next_free_slot;
                }
                if (next_free_slot < size_) {
                    used[next_free_slot] = true;
                    seeds_[bucket_idx] = -static_cast<std::int32_t>(next_free_slot + 1);
                }
                continue;
            }

            std::int32_t seed = 1;
            std::vector<std::size_t> chosen;
            chosen.reserve(b_size);

            while (true) {
                chosen.clear();
                bool collision = false;

                for (size_t i = b_start; i < b_end; ++i) {
                    const auto key_index = bucket_data[i];
                    const std::size_t slot = seededHash(keys[key_index], seed) % size_;
                    if (used[slot] || std::find(chosen.begin(), chosen.end(), slot) != chosen.end()) {
                        collision = true;
                        break;
                    }
                    chosen.push_back(slot);
                }

                if (!collision) {
                    for (const auto slot : chosen) {
                        used[slot] = true;
                    }
                    seeds_[bucket_idx] = seed;
                    break;
                }
                ++seed;
            }
        }
        std::cout << std::endl;
    }

    [[nodiscard]] std::size_t lookup(uint64_t key) const {
        if (size_ == 0) return 0;
        const std::size_t bucket = baseHash(key) % size_;
        g_cache_counter.record_miss(); // Access seeds_
        const std::int32_t seed = seeds_[bucket];
        if (seed < 0) {
            return static_cast<std::size_t>(-seed - 1);
        }
        return seededHash(key, seed) % size_;
    }

    [[nodiscard]] std::size_t size() const { return size_; }
    [[nodiscard]] std::size_t table_size() const { return size_; }

private:
    static std::uint64_t baseHash(uint64_t key) {
        // Standard 64-bit integer hash
        key = (key ^ (key >> 30)) * 0xbf58476d1ce4e5b9ULL;
        key = (key ^ (key >> 27)) * 0x94d049bb133111ebULL;
        key = key ^ (key >> 31);
        return key;
    }

    static std::uint64_t seededHash(uint64_t key, const std::int32_t seed) {
        return baseHash(key) ^ (0x9e3779b97f4a7c15ULL * static_cast<std::uint64_t>(seed));
    }

    std::size_t size_ = 0;
    std::vector<std::int32_t> seeds_;
};

} // namespace flat_index

#endif // MPHF_HPP
