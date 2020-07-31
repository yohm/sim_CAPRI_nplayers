# Code for the evaluation of CAPRI-n strategy

Here is the list of executables made by this projects.

- code to verify the three conditions using a graph-theoretic calculation
    - test_StrategyN2M3.out (for n=2)
    - test_StrategyN3M5.out (for n=3)
- evolutionary game simulation
    - main_evo_n2.out (for n=2)
    - main_evo_n3.out (for n=3)
- Monte Carlo test of defensibility against memory-1 strategies
    - main_run_N2_defensibility.out (for n=2)
    - main_run_N3_defensibility.out (for n=3)
    - main_run_N4_defensibility.out (for n=4)
- code for a unit test
    - test_DirectedGraph.out
    
# How to build

You need a C++ compiler supporting OpenMP and CMake.
The code is dependent on [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) library for linear algebraic calculation.

For a macOS user, you can prepare these prerequisites as follows.

```sh
brew install eigen
brew install libomp
brew install cmake
```

Then, create a build directory and run the cmake as follows.

```sh
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release
make
```

The executables are created in the "cpp/" directory.

