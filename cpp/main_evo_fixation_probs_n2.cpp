//
// Created by Yohsuke Murase on 2020/06/04.
//

#include <iostream>
#include <ostream>
#include <vector>
#include <array>
#include <random>
#include <cassert>
#include <omp.h>
#include <mpi.h>
#include <Eigen/Dense>
#include <fstream>
#include "StrategyN2M3.hpp"

// calculate the distribution of fixation probability rho
// against randomly selected N2M3 deterministic strategies

std::array<double,2> CalcPayoffs(const std::array<double,64>& stationary_state, double benefit) {
  const double cost = 1.0;
  std::array<double,2> ans = {0.0, 0.0};
  for (size_t i = 0; i < 64; i++) {
    StateN2M3 s(i);
    double pa = 0.0, pb = 0.0;
    if (s.a_1 == C) { pa -= cost; pb += benefit; }
    if (s.b_1 == C) { pb -= cost; pa += benefit; }
    ans[0] += stationary_state[i] * pa;
    ans[1] += stationary_state[i] * pb;
  }
  return std::move(ans);
}

double FixationProb(size_t N, double sigma, double e, double benefit, const StrategyN2M3 &res, const StrategyN2M3 &mut, double s_yy) {
  auto a_xx = mut.StationaryState(e);
  auto a_xy = mut.StationaryState(e, &res);

  double s_xx = CalcPayoffs(a_xx, benefit)[0];
  auto _xy = CalcPayoffs(a_xy, benefit);
  double s_xy = _xy[0];
  double s_yx = _xy[1];

  // \frac{1}{\rho} = \sum_{i=0}^{N-1} \exp\left( \sigma \sum_{j=1}^{i} \left[(N-j-1)s_{yy} + js_{yx} - (N-j)s_{xy} - (j-1)s_{xx} \right] \right) \\
  //                = \sum_{i=0}^{N-1} \exp\left( \frac{\sigma i}{2} \left[(-i+2N-3)s_{yy} + (i+1)s_{yx} - (-i+2N-1)s_{xy} - (i-1)s_{xx} \right] \right)

  double num_games = (N-1);
  s_xx /= num_games;
  s_yy /= num_games;
  s_xy /= num_games;
  s_yx /= num_games;
  double rho_inv = 0.0;
  for (int i=0; i < N; i++) {
    double x = sigma * i * 0.5 * (
        (2*N-3-i) * s_yy
            + (i+1) * s_yx
            - (2*N-1-i) * s_xy
            - (i-1) * s_xx
    );
    rho_inv += std::exp(x);
  }
  return 1.0 / rho_inv;
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

  std::uniform_int_distribution<uint64_t > dist(0, std::numeric_limits<uint64_t>::max() );

  const size_t NUM_BINS = 1000;
  std::vector<size_t> counts(NUM_BINS, 0ul);
  size_t robust_count = 0ul;

  int my_rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  int num_procs = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

  if (n_resident > 0) {
    #pragma omp parallel
    {
      std::vector<size_t> counts_tl(NUM_BINS, 0ul);
      size_t robust_count_tl = 0ul;
      int th = omp_get_thread_num();
      std::seed_seq seq = {seed, my_rank, th};
      std::mt19937_64 rnd_tl(seq);
      size_t my_n_resident = n_resident / num_procs;
      if (my_rank < n_resident % num_procs) { my_n_resident++; }
      // std::cerr << "my_n_resident: " << my_n_resident << " " << my_rank << ' ' << th << std::endl;

      #pragma omp for
      for (size_t i = 0; i < my_n_resident; i++) {
        uint64_t r = dist(rnd_tl);
        StrategyN2M3 res(r);
        auto a_yy = res.StationaryState(e);
        double s_yy = CalcPayoffs(a_yy, benefit)[0];
        for (size_t j = 0; j < n_mutants; j++) {
          uint64_t r2 = dist(rnd_tl);
          StrategyN2M3 mut(r2);
          double rho = FixationProb(N, sigma, e, benefit, res, mut, s_yy);
          size_t b = static_cast<size_t>(rho * NUM_BINS);
          counts_tl[b]++;
          if (rho <= 1.0 / N) { robust_count_tl++; }
        }
      }
      // copy from thread local variables
      for (size_t i = 0; i < counts_tl.size(); i++) {
        #pragma omp atomic update
        counts[i] += counts_tl[i];
      }
      #pragma omp atomic update
      robust_count += robust_count_tl;
    }
  }
  else {
    #pragma omp parallel
    {
      std::vector<size_t> counts_tl(NUM_BINS, 0ul);
      size_t robust_count_tl = 0ul;
      int th = omp_get_thread_num();
      std::seed_seq seq = {seed, my_rank, th};
      std::mt19937_64 rnd_tl(seq);

      StrategyN2M3 res = StrategyN2M3::CAPRI2();
      auto a_yy = res.StationaryState(e);
      double s_yy = CalcPayoffs(a_yy, benefit)[0];

      size_t my_n_mutants = n_mutants / num_procs;
      if (my_rank < n_mutants % num_procs) { my_n_mutants++; }
      //std::cerr << "my_n_mutants: " << my_n_mutants << ' ' << my_rank << ' ' << th << std::endl;

      #pragma omp for
      for (size_t j = 0; j < my_n_mutants; j++) {
        uint64_t r2 = dist(rnd_tl);
        StrategyN2M3 mut(r2);
        double rho = FixationProb(N, sigma, e, benefit, res, mut, s_yy);
        size_t b = static_cast<size_t>(rho * NUM_BINS);
        counts_tl[b] += 1;
        if (rho <= 1.0 / N) { robust_count_tl++; }
      }
      // copy from thread local variables
      for (size_t i = 0; i < counts_tl.size(); i++) {
      #pragma omp atomic update
        counts[i] += counts_tl[i];
      }
      #pragma omp atomic update
      robust_count += robust_count_tl;
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
