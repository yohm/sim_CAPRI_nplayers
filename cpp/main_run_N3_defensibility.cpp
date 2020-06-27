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
#include "StrategyN3M5.hpp"


class Mem1StrategyN3 {
 public:
  Mem1StrategyN3(double _c0, double _c1, double _c2, double _d0, double _d1, double _d2) :
      pc({_c0,_c1,_c2}), pd({_d0,_d1,_d2}) {};
  Mem1StrategyN3(std::mt19937_64 &rnd) {
    std::uniform_real_distribution<double> uni;
    for (size_t i = 0; i < 3; i++) {
      pc[i] = uni(rnd);
      pd[i] = uni(rnd);
    }
  }
  std::array<double,3> pc, pd;
  double CProb(const StateN3M5 &c) const {
    size_t nd = 0;
    if(c.hb[0]) nd++;
    if(c.hc[0]) nd++;
    return c.ha[0] ? pd[nd] : pc[nd];
  }
};

std::array<double,3> Run(size_t t_max, double benefit, std::mt19937_64 &rnd) {
  std::uniform_real_distribution<double> uni;

  StateN3M5 current(0);
  StrategyN3M5 strA = StrategyN3M5::CAPRI3();
  Mem1StrategyN3 strB(rnd);
  Mem1StrategyN3 strC(rnd);

  std::array<double,3> payoffs = {0.0, 0.0, 0.0};
  for (size_t t = 0; t < t_max; t++) {
    Action act_a = strA.ActionAt(current);
    double pb = strB.CProb(current.StateFromB());
    double pc = strC.CProb(current.StateFromC());
    Action act_b = (uni(rnd) < pb ? C : D);
    Action act_c = (uni(rnd) < pc ? C : D);
    current = current.NextState(act_a, act_b, act_c);

    std::array<double,3> p = {0.0, 0.0, 0.0};

    // public goods game with multiplication factor benefit
    // int n_c = 0;
    // if (act_a == C) {
    //   p[0] -= 1.0;
    //   n_c++;
    // }
    // if (act_b == C) {
    //   p[1] -= 1.0;
    //   n_c++;
    // }
    // if (act_c == C) {
    //   p[2] -= 1.0;
    //   n_c++;
    // }
    // for (int i = 0; i < 3; i++) { p[i] += benefit * n_c / 3.0; }

    // Donation game
    if (act_a == C) {
      p[0] -= 1.0;
      p[1] += benefit/2.0;
      p[2] += benefit/2.0;
    }
    if (act_b == C) {
      p[1] -= 1.0;
      p[0] += benefit/2.0;
      p[2] += benefit/2.0;
    }
    if (act_c == C) {
      p[2] -= 1.0;
      p[0] += benefit/2.0;
      p[1] += benefit/2.0;
    }
    double margin = 2.0 * (benefit/2.0 + 1.0);

    for (int i = 0; i < 3; i++) {
      payoffs[i] += p[i];
      if (i > 0 && payoffs[i] >= payoffs[0] + margin) {
        throw std::runtime_error("defensibility is violated");
      }
    }
  }

  for (int i = 0; i < 3; i++) { payoffs[i] /= t_max; }
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

  std::vector<std::array<double,3> > results(n_samples);

#pragma omp parallel for
  for (uint64_t n = 0; n < n_samples; n++) {
    if (n % 1000 == 0) { std::cerr << "n: " << n << std::endl; }
    std::seed_seq seed = {seed_base, n};
    std::mt19937_64 rnd(seed);
    results[n] = Run(t_max, benefit, rnd);
  }

  for (auto x: results) {
    std::cout << x[0] << ' ' << x[1] << ' ' << x[2] << std::endl;
  }
  return 0;
}

