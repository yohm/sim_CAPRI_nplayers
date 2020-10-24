//
// Created by Yohsuke Murase on 2020/06/04.
//

#include <iostream>
#include <ostream>
#include <vector>
#include <array>
#include <random>
#include <cassert>
#include <mpi.h>
#include <Eigen/Dense>
#include <fstream>
#include "StrategyN3M5.hpp"

// calculate the distribution of fixation probability rho
// against randomly selected N2M3 deterministic strategies

std::array<double,3> CalcPayoffs(const std::array<double,StrategyN3M5::N>& stationary_state, double benefit) {
  const double cost = 1.0;
  std::array<double,3> ans = {0.0, 0.0, 0.0};
  for (size_t i = 0; i < StrategyN3M5::N; i++) {
    StateN3M5 s(i);
    double pa = 0.0, pb = 0.0, pc = 0.0;
    if (s.ha[0] == false) { pa -= cost; pb += benefit / 2.0; pc += benefit / 2.0; }
    if (s.hb[0] == false) { pb -= cost; pa += benefit / 2.0; pc += benefit / 2.0; }
    if (s.hc[0] == false) { pc -= cost; pa += benefit / 2.0; pb += benefit / 2.0; }
    ans[0] += stationary_state[i] * pa;
    ans[1] += stationary_state[i] * pb;
    ans[2] += stationary_state[i] * pc;
  }
  return ans;
}

double FixationProb(size_t N, double sigma, double e, double benefit, const StrategyN3M5 &res, const StrategyN3M5 &mut, double s_yyy) {
  auto a_xxx = mut.StationaryState(e);
  auto a_xxy = mut.StationaryState(e, &mut, &res);
  auto a_xyy = mut.StationaryState(e, &res, &res);

  double s_xxx = CalcPayoffs(a_xxx, benefit)[0];
  auto _xxy = CalcPayoffs(a_xxy, benefit);
  double s_xxy = _xxy[0];
  double s_yxx = _xxy[2];
  auto _xyy = CalcPayoffs(a_xyy, benefit);
  double s_xyy = _xyy[0];
  double s_yyx = _xyy[1];

  // rho_inv = \sum_{i=0}^{N-1} exp(sigma[S]),
  // where S is defined as
  // S =  i/6(i^2−3iN+6i+3N^2−12N+11)s_{yyy}
  //     −i/6(i+1)(2i−3N+4)s_{yyx}
  //     +i/6(i−1)(i+1)s_{yxx}
  //     −i/6(i^2−3i(N−1)+3N^2−6N+2)s_{xyy}
  //     +i/6(i−1)(2i−3N+2))s_{xxy}
  //     −i/6(i^2−3i+2)s_{xxx}

  double num_games = (N-1) * (N-2) / 2.0;
  double rho_inv = 0.0;
  for (int i=0; i < N; i++) {
    double x = sigma / num_games * (i / 6.0) * (
        (i*i - 3.0*i*N + 6.0*i + 3.0*N*N - 12.0*N + 11.0) * s_yyy
            - (i+1.0) * (2.0*i - 3.0*N + 4) * s_yyx
            + (i-1.0) * (i+1.0) * s_yxx
            - (i*i - 3.0*i*(N-1.0) + 3.0*N*N - 6.0*N + 2.0) * s_xyy
            + (i-1.0) * (2.0*i-3.0*N+2.0) * s_xxy
            - (i*i - 3.0*i + 2.0) * s_xxx
    );
    rho_inv += std::exp(x);
  }

  return 1.0 / rho_inv;
}

std::bitset<StrategyN3M5::N> DrawRandomBit32768(std::mt19937_64 & rnd) {
  std::uniform_int_distribution<uint64_t > dist(0, std::numeric_limits<uint64_t>::max() );

  std::bitset<StrategyN3M5::N> ans;
  int rep = StrategyN3M5::N / 64;
  assert(StrategyN3M5::N % 64 == 0);
  for (int i = 0; i < rep; i++) {
    std::bitset<StrategyN3M5::N> mask = dist(rnd);
    mask <<= (64 * i);
    ans |= mask;
  }
  return ans;
}

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);
  Eigen::initParallel();
  if( argc != 8 ) {
    std::cerr << "Error : invalid argument" << std::endl;
    std::cerr << "  Usage : " << argv[0] << " <N> <sigma> <e> <benefit> <resident 0:caprin, +:num_resident_samples> <num_mutants> <seed>" << std::endl;
    MPI_Finalize();
    return 1;
  }

  size_t N = std::strtoul(argv[1], nullptr,0);
  double sigma = std::strtod(argv[2], nullptr);
  double e = std::strtod(argv[3], nullptr);
  double benefit = std::strtod(argv[4], nullptr);

  long n_resident = std::strtol(argv[5], nullptr, 0);
  long n_mutants = std::strtol(argv[6], nullptr, 0);
  int seed = std::strtol(argv[7], nullptr, 0);


  const size_t NUM_BINS = 1000;
  std::vector<size_t> counts(NUM_BINS, 0ul);
  size_t robust_count = 0ul;

  int my_rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  int num_procs = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

  if (n_resident > 0) {
    {
      std::seed_seq seq = {seed, my_rank};
      std::mt19937_64 rnd(seq);
      size_t my_n_resident = n_resident / num_procs;
      if (my_rank < n_resident % num_procs) { my_n_resident++; }
      // std::cerr << "my_n_resident: " << my_n_resident << " " << my_rank << ' ' << th << std::endl;

      for (size_t i = 0; i < my_n_resident; i++) {
        std::cerr << "i: " << i << " " << my_rank << std::endl;
        StrategyN3M5 res( DrawRandomBit32768(rnd) );
        auto a_yyy = res.StationaryState(e);
        double s_yyy = CalcPayoffs(a_yyy, benefit)[0];
        for (size_t j = 0; j < n_mutants; j++) {
          std::cerr << "j: " << j << " " << my_rank << std::endl;
          StrategyN3M5 mut( DrawRandomBit32768(rnd) );
          double rho = FixationProb(N, sigma, e, benefit, res, mut, s_yyy);
          std::cerr << "  rho: " << rho << std::endl;
          size_t b = static_cast<size_t>(rho * NUM_BINS);
          counts[b]++;
          if (rho <= 1.0 / N) { robust_count++; }
        }
      }
    }
  }
  else {
    {
      std::seed_seq seq = {seed, my_rank};
      std::mt19937_64 rnd(seq);

      StrategyN3M5 res = StrategyN3M5::CAPRI3();
      auto a_yyy = res.StationaryState(e);
      double s_yyy = CalcPayoffs(a_yyy, benefit)[0];

      size_t my_n_mutants = n_mutants / num_procs;
      if (my_rank < n_mutants % num_procs) { my_n_mutants++; }
      //std::cerr << "my_n_mutants: " << my_n_mutants << ' ' << my_rank << ' ' << th << std::endl;

      for (size_t j = 0; j < my_n_mutants; j++) {
        StrategyN3M5 mut( DrawRandomBit32768(rnd) );
        double rho = FixationProb(N, sigma, e, benefit, res, mut, s_yyy);
        size_t b = static_cast<size_t>(rho * NUM_BINS);
        counts[b] += 1;
        if (rho <= 1.0 / N) { robust_count++; }
      }
    }
  }

  // reduce counts
  std::vector<size_t> all_counts(NUM_BINS, 0ul);
  size_t all_robust_count = 0ul;
  MPI_Reduce(counts.data(), all_counts.data(), all_counts.size(), MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&robust_count, &all_robust_count, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

  // print counts
  if (my_rank == 0) {
    std::ofstream fout("dist.dat");
    double dx = 1.0 / NUM_BINS;
    double total = (n_resident == 0) ? n_mutants : n_resident * n_mutants;
    for (size_t i = 0; i < NUM_BINS; i++) {
      fout << i * dx << ' ' << (double)all_counts[i]/total << std::endl;
    }

    std::cerr << "robust_count/total: " << all_robust_count << " / " << total << " : " << (double)all_robust_count/total << std::endl;
  }

  MPI_Finalize();
  return 0;
}
