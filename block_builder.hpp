#ifndef BLOCK_BUILDER_HPP
#define BLOCK_BUILDER_HPP

#include <vector>
#include <cstddef>
#include <cstdint>
#include "data_model.hpp"
#include "flat_block.hpp"

namespace flat_index {

struct PackedBlocks {
    std::vector<std::byte> bytes;            // contiguous, multiple of 64
    // Note: first_block_byte_offset is NOT known here, it's known by the assembler
};

PackedBlocks pack_unitig(const InUnitig& u,
                         const std::vector<InOccurrence>& utab_rows,
                         const std::vector<InEdge>& etab_rows,
                         int k);

} // namespace flat_index

#endif // BLOCK_BUILDER_HPP
