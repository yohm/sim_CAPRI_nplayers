///
// Created by Yohsuke Murase on 2020/06/24.
//

#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <bitset>
#include <random>
#include <cstdint>
#include "Action.hpp"
#include "StrategyN2M3.hpp"


class Mem1StrategyN2 {
 public:
  Mem1StrategyN2(double _cc, double _cd, double _dc, double _dd) : probs({_cc, _cd, _dc, _dd}) {};
  Mem1StrategyN2(std::mt19937_64 &rnd) {
    std::uniform_real_distribution<double> uni;
    for (size_t i = 0; i < 4; i++) { probs[i] = uni(rnd); }
  }
  std::array<double,4> probs;
  double CProb(const StateN2M3 &c) const {
    if (c.a_1 == C) {
      if (c.b_1 == C) { return probs[0]; }
      else { return probs[1]; }
    }
    else {
      if (c.b_1 == C) { return probs[2]; }
      else { return probs[3]; }
    }
  }
};

std::array<double,2> Run(size_t t_max, double benefit, std::mt19937_64 &rnd) {
  std::uniform_real_distribution<double> uni;

  StateN2M3 current(0ull);
  // StrategyN2M3 strA = StrategyN2M3::CAPRI();
  // StrategyN2M3 strA = StrategyN2M3::sCAPRI2();
  StrategyN2M3 strA = StrategyN2M3::CAPRI2();
  Mem1StrategyN2 strB(rnd);

  std::array<double,2> payoffs = {0.0, 0.0};
  for (size_t t = 0; t < t_max; t++) {
    Action act_a = strA.ActionAt(current);
    double pb = strB.CProb(current.SwapAB());
    Action act_b = (uni(rnd) < pb ? C : D);
    current = current.NextState(act_a, act_b);

    std::array<double,2> p = {0.0, 0.0};
    int n_c = 0;
    if (act_a == C) {
      p[0] -= 1.0;
      n_c++;
    }
    if (act_b == C) {
      p[1] -= 1.0;
      n_c++;
    }
    for (int i = 0; i < 2; i++) { p[i] += benefit * n_c / 2.0; }

    for (int i = 0; i < 2; i++) {
      payoffs[i] += p[i];
      if (i > 0 && payoffs[i] >= payoffs[0] + 2.0) {
        throw std::runtime_error("defensibility is violated");
      }
    }
  }

  for (int i = 0; i < 2; i++) { payoffs[i] /= t_max; }
  return payoffs;
}

int main(int argc, char* argv[]) {
  if (argc != 5) {
    std::cerr << "invalid number of arguments" << std::endl;
    std::cerr << "usage: " << argv[0] << " <t_max> <benefit> <n_samples> <seed>" << std::endl;
  }

  uint64_t t_max = std::strtoull(argv[1], nullptr,0);
  double benefit = std::strtod(argv[2], nullptr);
  uint64_t n_samples = std::strtoull(argv[3], nullptr,0);
  uint64_t seed_base = std::strtoull(argv[4], nullptr,0);

  std::vector<std::array<double,2> > results(n_samples);

#pragma omp parallel for
  for (uint64_t n = 0; n < n_samples; n++) {
    if (n % 1000 == 0) { std::cerr << "n: " << n << std::endl; }
    std::seed_seq seed = {seed_base, n};
    std::mt19937_64 rnd(seed);
    results[n] = Run(t_max, benefit, rnd);
  }

  for (auto x: results) {
    std::cout << x[0] << ' ' << x[1] << std::endl;
  }
  return 0;
}

