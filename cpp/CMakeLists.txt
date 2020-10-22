cmake_minimum_required(VERSION 3.7)
project(cpp)

set(CMAKE_CXX_STANDARD 14)

#set(CMAKE_VERBOSE_MAKEFILE 1)

set(SOURCE_N2M3_FILES Action.cpp Action.hpp StrategyN2M3.cpp StrategyN2M3.hpp DirectedGraph.hpp DirectedGraph.cpp UnionFind.hpp Partition.hpp)
include_directories(SYSTEM /usr/local/include/eigen3)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(test_StrategyN2M3.out test_StrategyN2M3.cpp ${SOURCE_N2M3_FILES})
add_executable(test_DirectedGraph.out test_DirectedGraph.cpp DirectedGraph.hpp DirectedGraph.cpp)

set(SOURCE_N3M5_FILES Action.cpp Action.hpp StrategyN3M5.cpp StrategyN3M5.hpp DirectedGraph.hpp DirectedGraph.cpp UnionFind.hpp Partition.hpp)
add_executable(test_StrategyN3M5.out test_StrategyN3M5.cpp ${SOURCE_N3M5_FILES})

add_executable(main_evo_n2.out main_evo_n2.cpp ${SOURCE_N2M3_FILES})
add_executable(main_evo_fixation_probs_n2.out main_evo_fixation_probs_n2.cpp ${SOURCE_N2M3_FILES})

find_package(OpenMP REQUIRED)
add_executable(main_evo_n3.out main_evo_n3.cpp ${SOURCE_N3M5_FILES})
target_link_libraries(main_evo_n3.out PRIVATE OpenMP::OpenMP_CXX)

add_executable(main_run_N4_defensibility.out main_run_N4_defensibility.cpp Action.cpp Action.hpp)
target_link_libraries(main_run_N4_defensibility.out PRIVATE OpenMP::OpenMP_CXX)
add_executable(main_run_N3_defensibility.out main_run_N3_defensibility.cpp ${SOURCE_N3M5_FILES})
target_link_libraries(main_run_N3_defensibility.out PRIVATE OpenMP::OpenMP_CXX)
add_executable(main_run_N2_defensibility.out main_run_N2_defensibility.cpp ${SOURCE_N2M3_FILES})
target_link_libraries(main_run_N2_defensibility.out PRIVATE OpenMP::OpenMP_CXX)

add_executable(main_run_N4_efficiency.out main_run_N4_efficiency.cpp Action.cpp Action.hpp)
target_link_libraries(main_run_N4_efficiency.out PRIVATE OpenMP::OpenMP_CXX)
