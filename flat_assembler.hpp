#ifndef FLAT_ASSEMBLER_HPP
#define FLAT_ASSEMBLER_HPP

#include <vector>
#include <unordered_map>
#include <cstddef>
#include <memory>
#include <cstdlib>
#include "data_model.hpp"

namespace flat_index {

template <typename T, std::size_t Alignment>
struct AlignedAllocator {
    using value_type = T;
    AlignedAllocator() noexcept = default;
    template <typename U> AlignedAllocator(const AlignedAllocator<U, Alignment>&) noexcept {}
    template <typename U> struct rebind { using other = AlignedAllocator<U, Alignment>; };
    T* allocate(std::size_t n) {
        if (n == 0) return nullptr;
        void* ptr = nullptr;
#if defined(_MSC_VER) || defined(__MINGW32__)
        ptr = _aligned_malloc(n * sizeof(T), Alignment);
        if (ptr == nullptr) throw std::bad_alloc();
#else
        if (posix_memalign(&ptr, Alignment, n * sizeof(T)) != 0) throw std::bad_alloc();
#endif
        return static_cast<T*>(ptr);
    }
    void deallocate(T* p, std::size_t) noexcept {
#if defined(_MSC_VER) || defined(__MINGW32__)
        _aligned_free(p);
#else
        free(p);
#endif
    }
    bool operator==(const AlignedAllocator&) const { return true; }
    bool operator!=(const AlignedAllocator&) const { return false; }
};

struct FlatIndex {
    std::vector<std::byte, AlignedAllocator<std::byte, 64>> flat;
    int k;
};

FlatIndex assemble(const std::vector<InUnitig>& unitigs,
                   const std::vector<InOccurrence>& occs,
                   const std::vector<InEdge>& edges,
                   int k);

} // namespace flat_index

#endif // FLAT_ASSEMBLER_HPP
