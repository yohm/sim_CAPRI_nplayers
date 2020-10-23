//
// Created by Yohsuke Murase on 2020/06/04.
//

#include <iostream>
#include <ostream>
#include <vector>
#include <array>
#include <string>
#include <map>
#include <random>
#include <cassert>
#include <fstream>
#include <Eigen/Dense>
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
  Eigen::initParallel();
  if( argc != 8 ) {
    std::cerr << "Error : invalid argument" << std::endl;
    std::cerr << "  Usage : " << argv[0] << " <N> <sigma> <e> <benefit> <resident 0:caprin, +:num_resident_samples> <num_mutants> <seed>" << std::endl;
    return 1;
  }

  size_t N = std::strtoul(argv[1], nullptr,0);
  double sigma = std::strtod(argv[2], nullptr);
  double e = std::strtod(argv[3], nullptr);
  double benefit = std::strtod(argv[4], nullptr);

  long n_resident = std::strtol(argv[5], nullptr, 0);
  long n_mutants = std::strtol(argv[6], nullptr, 0);
  uint64_t seed = std::strtoull(argv[7], nullptr, 0);

  std::mt19937_64 rnd(seed);
  std::uniform_int_distribution<uint64_t > dist(0, std::numeric_limits<uint64_t>::max() );

  const size_t NUM_BINS = 100;
  std::vector<size_t> counts(NUM_BINS, 0ul);
  size_t robust_count = 0ul;

  if (n_resident > 0) {
    for (size_t i = 0; i < n_resident; i++) {
      uint64_t r = dist(rnd);
      StrategyN2M3 res(r);
      auto a_yy = res.StationaryState(e);
      double s_yy = CalcPayoffs(a_yy, benefit)[0];
      for (size_t j = 0; j < n_mutants; j++) {
        uint64_t r2 = dist(rnd);
        StrategyN2M3 mut(r2);
        double rho = FixationProb(N, sigma, e, benefit, res, mut, s_yy);
        size_t b = static_cast<size_t>(rho * NUM_BINS);
        counts[b] += 1;
        if (rho <= 1.0 / N) { robust_count++; }
      }
    }
  }
  else {
    StrategyN2M3 res = StrategyN2M3::CAPRI2();
    auto a_yy = res.StationaryState(e);
    double s_yy = CalcPayoffs(a_yy, benefit)[0];
    for (size_t j = 0; j < n_mutants; j++) {
      uint64_t r2 = dist(rnd);
      StrategyN2M3 mut(r2);
      double rho = FixationProb(N, sigma, e, benefit, res, mut, s_yy);
      size_t b = static_cast<size_t>(rho * NUM_BINS);
      counts[b] += 1;
      if (rho <= 1.0 / N) { robust_count++; }
    }
  }

  // print counts
  double dx = 1.0 / NUM_BINS;
  double total = (n_resident == 0) ? n_mutants : n_resident * n_mutants;
  for (size_t i = 0; i < NUM_BINS; i++) {
    std::cout << i * dx << ' ' << (double)counts[i]/total << std::endl;
  }

  std::cerr << "robust_count/total: " << robust_count << " / " << total << " : " << (double)robust_count/total << std::endl;

  return 0;
}
