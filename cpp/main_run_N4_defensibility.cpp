//
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

class StateN4M7 {
 public:
  StateN4M7() : ha(0), hb(0), hc(0), hd(0) {}
  std::bitset<7> ha,hb,hc,hd;
  void Update(Action a, Action b, Action c, Action d) {
    ha <<= 1; hb <<= 1; hc <<= 1; hd <<= 1;
    if (a == D) { ha.set(0); }
    if (b == D) { hb.set(0); }
    if (c == D) { hc.set(0); }
    if (d == D) { hd.set(0); }
  }
  std::array<double,4> Payoffs(double b) const {
    std::array<double,4> ans = {0.0, 0.0, 0.0, 0.0};
    // Public goods game with multiplication factor b
    // int n_c = 0;
    // if (ha[0] == C) { n_c++; ans[0] = -1.0; }
    // if (hb[0] == C) { n_c++; ans[1] = -1.0; }
    // if (hc[0] == C) { n_c++; ans[2] = -1.0; }
    // if (hd[0] == C) { n_c++; ans[3] = -1.0; }
    // for (size_t i = 0; i < 4; i++) {
    //   ans[i] += n_c * b / 4.0;
    // }

    // Donation game with b
    if (!ha[0]) { ans[0] -= 1.0; ans[1] += b / 3.0; ans[2] += b / 3.0; ans[3] += b / 3.0; }
    if (!hb[0]) { ans[1] -= 1.0; ans[2] += b / 3.0; ans[3] += b / 3.0; ans[0] += b / 3.0; }
    if (!hc[0]) { ans[2] -= 1.0; ans[3] += b / 3.0; ans[0] += b / 3.0; ans[1] += b / 3.0; }
    if (!hd[0]) { ans[3] -= 1.0; ans[0] += b / 3.0; ans[1] += b / 3.0; ans[2] += b / 3.0; }

    return std::move(ans);
  }
};

Action CAPRI4_Action(StateN4M7 c) {
  {  // when Nd >= 4, prescribe D
    size_t na = c.ha.count();
    size_t nb = c.hb.count();
    size_t nc = c.hc.count();
    size_t nd = c.hd.count();
    size_t l = std::max(na, std::max(nb, std::max(nc, nd)));
    size_t s = std::min(na, std::min(nb, std::min(nc, nd)));
    if (l - s >= 4) { return D; }
  }

  // find last_cccc
  size_t last_cccc = 7;
  {
    std::bitset<7> h_or = c.ha | c.hb | c.hc | c.hd;
    for (size_t t = 0; t < 7; t++) {
      if (h_or[t] == false) {  // CCCC is found at t step before
        last_cccc = t;
        break;
      }
    }
  }

  // C: last action is mutual cooperation
  if (last_cccc == 0) {
    return C;
  }
  else if (last_cccc > 0 && last_cccc < 7) {
    std::bitset<7> mask = 0;
    for (size_t t = 0; t < last_cccc; t++) { mask.set(t); }
    // A: Accept punishment by prescribing *C* if all your relative payoffs are at least zero.
    size_t pa = (c.ha & mask).count();
    size_t pb = (c.hb & mask).count();
    size_t pc = (c.hc & mask).count();
    size_t pd = (c.hd & mask).count();
    if (pa >= pb && pa >= pc && pa >= pd) {
      return C;
    }
    // P: Punish by *D* if any of your relative payoffs is negative.
    else {
      return D;
    }
  }

  // R: grab the chance to recover
  // 11...0, 11...0, 11...0, 11...1 or its permutations
  {
    std::bitset<7> h_and = c.ha & c.hb & c.hc & c.hd;
    if (h_and.count() == 6 && h_and[0] == false) {
      size_t n = 0;
      if (c.ha[0]) n++;
      if (c.hb[0]) n++;
      if (c.hc[0]) n++;
      if (c.hd[0]) n++;
      if (n == 1) return C;
    }
  }

  return D;
}


class Mem1Strategy {
 public:
  Mem1Strategy(double _c0, double _c1, double _c2, double _c3, double _d0, double _d1, double _d2, double _d3) :
  pc({_c0,_c1,_c2,_c3}), pd({_d0,_d1,_d2,_d3}) {};
  Mem1Strategy(std::mt19937_64 &rnd) {
    std::uniform_real_distribution<double> uni;
    for (size_t i = 0; i < 4; i++) {
      pc[i] = uni(rnd);
      pd[i] = uni(rnd);
    }
  }
  std::array<double,4> pc, pd;
  double CProb(const StateN4M7 &c, int x) const { // x=0 => B, x=1 => C, x=2 => D
    if (x == 0) {
      size_t nd = 0;
      if(c.ha[0]) nd++;
      if(c.hc[0]) nd++;
      if(c.hd[0]) nd++;
      return c.hb[0] ? pd[nd] : pc[nd];
    }
    else if (x == 1) {
      size_t nd = 0;
      if(c.ha[0]) nd++;
      if(c.hb[0]) nd++;
      if(c.hd[0]) nd++;
      return c.hc[0] ? pd[nd] : pc[nd];
    }
    else {
      size_t nd = 0;
      if(c.ha[0]) nd++;
      if(c.hb[0]) nd++;
      if(c.hc[0]) nd++;
      return c.hd[0] ? pd[nd] : pc[nd];
    }
  }
};

std::array<double,4> Run(size_t t_max, double benefit, std::mt19937_64 &rnd) {
  std::uniform_real_distribution<double> uni;

  StateN4M7 current;
  Mem1Strategy strB(rnd);
  Mem1Strategy strC(rnd);
  Mem1Strategy strD(rnd);

  std::array<double,4> payoffs = {0.0, 0.0, 0.0, 0.0};
  for (size_t t = 0; t < t_max; t++) {
    Action act_a = CAPRI4_Action(current);
    double pb = strB.CProb(current, 0);
    double pc = strC.CProb(current, 1);
    double pd = strD.CProb(current, 2);
    Action act_b = (uni(rnd) < pb ? C : D);
    Action act_c = (uni(rnd) < pc ? C : D);
    Action act_d = (uni(rnd) < pd ? C : D);
    StateN4M7 old = current;
    current.Update(act_a, act_b, act_c, act_d);
    auto p = current.Payoffs(benefit);
    for (int i = 0; i < 4; i++) {
      payoffs[i] += p[i];
      if (i > 0 && payoffs[i] >= payoffs[0] + 2.0*(benefit+1.0) ) {
        throw std::runtime_error("defensibility is violated");
      }
    }
  }

  for (int i = 0; i < 4; i++) { payoffs[i] /= t_max; }
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

  std::vector<std::array<double,4> > ans(n_samples);
#pragma omp parallel for
  for (uint64_t n = 0; n < n_samples; n++) {
    if (n % 1000 == 0) { std::cerr << "n: " << n << std::endl; }
    std::seed_seq seed = {seed_base, n};
    std::mt19937_64 rnd(seed);
    ans[n] = Run(t_max, benefit, rnd);
  }

  for (auto x: ans) {
    std::cout << x[0] << ' ' << x[1] << ' ' << x[2] << ' ' << x[3] << std::endl;
  }
  return 0;
}
