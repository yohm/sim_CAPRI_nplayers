cmake_minimum_required(VERSION 3.7)
project(cpp)

set(CMAKE_CXX_STANDARD 11)

#set(CMAKE_VERBOSE_MAKEFILE 1)

set(SOURCE_FILES Action.cpp Action.hpp StrategyN2M3.cpp StrategyN2M3.hpp DirectedGraph.hpp DirectedGraph.cpp UnionFind.hpp)
include_directories(SYSTEM /usr/local/include/eigen3)

add_executable(test_StrategyN2M3 test_StrategyN2M3.cpp ${SOURCE_FILES})
add_executable(test_DirectedGraph test_DirectedGraph.cpp DirectedGraph.hpp DirectedGraph.cpp)

