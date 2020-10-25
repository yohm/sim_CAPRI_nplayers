#!/bin/bash -eux

#mpiFCCpx -Nclang -std=c++14 -Kfast -o main_evo_fixation_probs_n3.out Action.cpp DirectedGraph.cpp StrategyN3M5.cpp main_evo_fixation_probs_n3.cpp -I $HOME/sandbox/eigen-3.3.7
mpiFCCpx -Nclang -std=c++14 -Kfast -o main_evo_fixation_probs_n2.out Action.cpp DirectedGraph.cpp StrategyN2M3.cpp main_evo_fixation_probs_n2.cpp -I $HOME/sandbox/eigen-3.3.7
#mpiFCCpx -Nclang -std=c++14 -Kfast -Kopenmp -o main_evo_fixation_probs_n3.out Action.cpp DirectedGraph.cpp StrategyN3M5.cpp main_evo_fixation_probs_n3.cpp -I $HOME/sandbox/eigen-3.3.7
#mpiFCCpx -Nclang -std=c++14 -Kfast -Kopenmp -o main_evo_fixation_probs_n3.out Action.cpp DirectedGraph.cpp StrategyN3M5.cpp main_evo_fixation_probs_n3.cpp -I $HOME/sandbox/eigen-3.3.7 -D EIGEN_NO_ASSERTION_CHECKING

# to avoid a bug in Eigen reported at https://bugs.archlinux.org/task/68118
# - define 'EIGEN_NO_ASSERTION_CHECKING' macro
# - comment out eigen-3.3.8/Eigen/src/Core/products/Parallelizer.h:162:40
#    `if (errorCount) EIGEN_THROW_X(Eigen::eigen_assert_exception());`

