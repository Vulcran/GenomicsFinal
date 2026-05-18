#ifndef SERIALIZER_HPP
#define SERIALIZER_HPP

#include "query.hpp"
#include <string>

namespace flat_index {

class IndexSerializer {
public:
    // Write a built FullIndex to a binary file.
    static void save(const FullIndex& idx, const std::string& path);

    // Read a FullIndex back from a binary file saved by save().
    static FullIndex load(const std::string& path);
};

} // namespace flat_index

#endif // SERIALIZER_HPP
