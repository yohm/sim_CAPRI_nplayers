cmake_minimum_required(VERSION 3.7)
project(cpp)

set(CMAKE_CXX_STANDARD 11)

#set(CMAKE_VERBOSE_MAKEFILE 1)

set(SOURCE_FILES Action.cpp Action.hpp Strategy.cpp Strategy.hpp DirectedGraph.hpp DirectedGraph.cpp)
include_directories(SYSTEM /usr/local/include/eigen3)

add_executable(test_Strategy test_Strategy.cpp ${SOURCE_FILES})
add_executable(test_DirectedGraph test_DirectedGraph.cpp ${SOURCE_FILES})

