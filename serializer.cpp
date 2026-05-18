#include "serializer.hpp"
#include <fstream>
#include <stdexcept>
#include <cstdint>

namespace flat_index {

// Binary format:
//   [4B magic] [4B version] [4B k]
//   [8B flat size] [flat bytes]
//   [8B pos_table count] [pos_table uint64s]
//   [MPHF: 8B size_ | 8B seeds count | int32_t seeds...]

static constexpr uint32_t MAGIC   = 0x464C4944; // "FLID"
static constexpr uint32_t VERSION = 1;

void IndexSerializer::save(const FullIndex& idx, const std::string& path) {
    std::ofstream out(path, std::ios::binary);
    if (!out) throw std::runtime_error("cannot write: " + path);

    auto w32 = [&](uint32_t v) { out.write(reinterpret_cast<const char*>(&v), 4); };
    auto w64 = [&](uint64_t v) { out.write(reinterpret_cast<const char*>(&v), 8); };

    w32(MAGIC);
    w32(VERSION);
    w32(static_cast<uint32_t>(idx.flat_idx.k));

    // Flat array.
    uint64_t flat_size = idx.flat_idx.flat.size();
    w64(flat_size);
    out.write(reinterpret_cast<const char*>(idx.flat_idx.flat.data()),
              static_cast<std::streamsize>(flat_size));

    // POS table.
    uint64_t pos_count = idx.pos_table.size();
    w64(pos_count);
    out.write(reinterpret_cast<const char*>(idx.pos_table.data()),
              static_cast<std::streamsize>(pos_count * 8));

    // MPHF.
    idx.mphf.save(out);

    if (!out) throw std::runtime_error("write error: " + path);
}

FullIndex IndexSerializer::load(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) throw std::runtime_error("cannot read: " + path);

    auto r32 = [&]() {
        uint32_t v;
        in.read(reinterpret_cast<char*>(&v), 4);
        return v;
    };
    auto r64 = [&]() {
        uint64_t v;
        in.read(reinterpret_cast<char*>(&v), 8);
        return v;
    };

    if (r32() != MAGIC)   throw std::runtime_error("bad magic in: " + path);
    if (r32() != VERSION) throw std::runtime_error("unsupported version in: " + path);

    FullIndex idx;
    idx.flat_idx.k = static_cast<int>(r32());

    // Flat array.
    uint64_t flat_size = r64();
    idx.flat_idx.flat.resize(flat_size);
    in.read(reinterpret_cast<char*>(idx.flat_idx.flat.data()),
            static_cast<std::streamsize>(flat_size));

    // POS table.
    uint64_t pos_count = r64();
    idx.pos_table.resize(pos_count);
    in.read(reinterpret_cast<char*>(idx.pos_table.data()),
            static_cast<std::streamsize>(pos_count * 8));

    // MPHF.
    idx.mphf.load(in);

    if (!in) throw std::runtime_error("read error or truncated file: " + path);
    return idx;
}

} // namespace flat_index
