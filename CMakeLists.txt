cmake_minimum_required(VERSION 3.15)
project(GenomicsFinal)

set(CMAKE_CXX_STANDARD 20)

# Core library sources shared by all targets.
set(CORE_SRCS
    block_builder.cpp
    flat_assembler.cpp
    query.cpp
    metrics.cpp
    loader.cpp
    serializer.cpp
)

add_executable(validate
    validate.cpp
    ${CORE_SRCS}
)

add_executable(bench
    bench.cpp
    ${CORE_SRCS}
)

add_executable(build_unitigs build_unitigs.cpp)
add_executable(build_occ     build_occ.cpp)

# Step 11/12: build a flat index from TSV or GFA input.
add_executable(flatindex_build
    flatindex_build.cpp
    ${CORE_SRCS}
)

# Step 11/12: query a serialized flat index.
add_executable(flatindex_query
    flatindex_query.cpp
    ${CORE_SRCS}
)
