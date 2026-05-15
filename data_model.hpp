#ifndef DATA_MODEL_HPP
#define DATA_MODEL_HPP

#include <cstdint>
#include <string>
#include <vector>

namespace flat_index {

struct InUnitig {
    uint32_t id;
    std::string seq;
};

struct InOccurrence {
    uint32_t unitig_id;
    uint16_t ref_id;
    uint32_t ref_pos;
    uint8_t entry_off;
    uint8_t walk_len;
    char orient; // '+' or '-'
};

struct InEdge {
    uint32_t from_id;
    char from_orient; // '+' or '-'
    uint32_t to_id;
    char to_orient; // '+' or '-'
    uint8_t overlap;
};

} // namespace flat_index

#endif // DATA_MODEL_HPP
