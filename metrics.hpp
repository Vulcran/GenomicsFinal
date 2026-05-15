#ifndef METRICS_HPP
#define METRICS_HPP

#include <cstddef>

namespace flat_index {

struct CacheCounter {
    size_t misses = 0;
    
    void record_miss() {
        misses++;
    }
    
    void reset() {
        misses = 0;
    }
};

extern CacheCounter g_cache_counter;

} // namespace flat_index

#endif // METRICS_HPP
