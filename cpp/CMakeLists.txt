cmake_minimum_required(VERSION 3.7)
project(cpp)

set(CMAKE_CXX_STANDARD 11)

#set(CMAKE_VERBOSE_MAKEFILE 1)

set(SOURCE_FILES Action.cpp Action.hpp Strategy.cpp Strategy.hpp DirectedGraph.hpp DirectedGraph.cpp TraceGSNegatives.hpp TraceGSNegatives.cpp CheckEfficiencyByGraph.hpp CheckEfficiencyByGraph.cpp)
find_package(MPI)
include_directories(SYSTEM ${MPI_INCLUDE_PATH} /usr/local/include/eigen3)

add_executable(test_Strategy test_Strategy.cpp ${SOURCE_FILES})
add_executable(test_DirectedGraph test_DirectedGraph.cpp ${SOURCE_FILES})

add_executable(main_trace_gSS_from_cccccd main_trace_gSS_from_cccccd.cpp ${SOURCE_FILES})
target_link_libraries(main_trace_gSS_from_cccccd ${MPI_LIBRARIES})

add_executable(main_trace_gS_from_negatives main_trace_gS_from_negatives.cpp ${SOURCE_FILES})
target_link_libraries(main_trace_gS_from_negatives ${MPI_LIBRARIES})

add_executable(main_filter_efficient_defensible main_filter_efficient_defensible.cpp ${SOURCE_FILES})
target_link_libraries(main_filter_efficient_defensible ${MPI_LIBRARIES})

add_executable(main_filter_distinguishable main_filter_distinguishable.cpp ${SOURCE_FILES})
target_link_libraries(main_filter_distinguishable ${MPI_LIBRARIES})

add_executable(main_check_one_by_one main_check_one_by_one.cpp ${SOURCE_FILES})
target_link_libraries(main_check_one_by_one ${MPI_LIBRARIES})

add_executable(main_print_ITG main_print_ITG.cpp ${SOURCE_FILES})
target_link_libraries(main_print_ITG ${MPI_LIBRARIES})

add_executable(count_num_strategies count_num_strategies.cpp)

add_executable(main_evolutionary_game main_evolutionary_game.cpp)
add_executable(main_payoff_histo main_payoff_histo.cpp)

