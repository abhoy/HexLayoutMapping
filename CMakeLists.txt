cmake_minimum_required(VERSION 3.10)
project(HexMapping 
    VERSION 1.0
    DESCRIPTION "Mapping Quantum Circuit to a Hexagonal Architecture"
    LANGUAGES CXX)


add_definitions(-std=c++11)
set(CORE_SOURCES 
    src/core/gate.hpp 
    src/core/gate.cpp
    src/core/circuit.hpp
    src/core/circuit.cpp
    src/core/graph.hpp
    src/core/graph.cpp)

set(ALGO_SOURCES 
    src/lib/global.hpp 
    src/lib/global.cpp 
    src/lib/local.hpp 
    src/lib/local.cpp 
    src/lib/optimize.hpp 
    src/lib/optimize.cpp)



set(TARGET_HEX_QX_MAPPER 
    src/Hx_mapper.cpp)

set(TARGET_HEX_QX_MAPPER_II 
    src/Hx_mapper_II.cpp)

add_library(Core SHARED ${CORE_SOURCES})

add_library(Algo SHARED ${ALGO_SOURCES})

target_include_directories(Core PUBLIC src/core)
target_include_directories(Algo PUBLIC src/lib)


add_executable(HEX_QxMapping    
    ${TARGET_HEX_QX_MAPPER}
    )
target_link_libraries(HEX_QxMapping Core Algo)


add_executable(HEX_QxMapping_II    
    ${TARGET_HEX_QX_MAPPER_II}
    )
target_link_libraries(HEX_QxMapping_II Core Algo)
