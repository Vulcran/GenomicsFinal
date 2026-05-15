cmake_minimum_required(VERSION 4.2)
project(GenomicsFinal)

set(CMAKE_CXX_STANDARD 20)

add_executable(GenomicsFinal main.cpp)

add_executable(validate 
    validate.cpp 
    block_builder.cpp 
    flat_assembler.cpp 
    query.cpp
    metrics.cpp
)

add_executable(bench 
    bench.cpp 
    block_builder.cpp 
    flat_assembler.cpp 
    query.cpp
    metrics.cpp
)
