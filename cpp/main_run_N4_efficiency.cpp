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
#include <utility>
#include "Action.hpp"

class StateN4M7 {
 public:
  StateN4M7(uint64_t i = 0) : ha(i), hb(i), hc(i), hd(i) {};
  std::bitset<7> ha,hb,hc,hd;
  void Update(Action a, Action b, Action c, Action d) {
    ha <<= 1; hb <<= 1; hc <<= 1; hd <<= 1;
    if (a == D) { ha.set(0); }
    if (b == D) { hb.set(0); }
    if (c == D) { hc.set(0); }
    if (d == D) { hd.set(0); }
  }
  bool Any() const {
    return (ha.any() || hb.any() || hc.any() || hd.any());
  }
  bool All() const {
    return (ha.all() && hb.all() && hc.all() && hd.all());
  }
  bool State1() const {
    unsigned long x = ha.to_ulong() + hb.to_ulong() + hc.to_ulong() + hd.to_ulong();
    return (x == 1);
  }
  std::array<double,4> Payoffs(double b) const {
    int n_c = 0;
    std::array<double,4> ans = {0.0, 0.0, 0.0, 0.0};
    if (ha[0] == C) { n_c++; ans[0] = -1.0; }
    if (hb[0] == C) { n_c++; ans[1] = -1.0; }
    if (hc[0] == C) { n_c++; ans[2] = -1.0; }
    if (hd[0] == C) { n_c++; ans[3] = -1.0; }
    for (size_t i = 0; i < 4; i++) {
      ans[i] += n_c * b / 4.0;
    }
    return std::move(ans);
  }
};

Action CAPRI4_Action(StateN4M7 c, int i) {
  std::bitset<7> ha, hb, hc, hd;
  if (i == 0) { ha = c.ha; hb = c.hb; hc = c.hc; hd = c.hd; }
  else if (i == 1) { ha = c.hb; hb = c.hc; hc = c.hd; hd = c.ha; }
  else if (i == 2) { ha = c.hc; hb = c.hd; hc = c.ha; hd = c.hb; }
  else { ha = c.hd; hb = c.ha; hc = c.hb; hd = c.hc; }

  {  // when Nd >= 4, prescribe D
    size_t na = ha.count();
    size_t nb = hb.count();
    size_t nc = hc.count();
    size_t nd = hd.count();
    size_t l = std::max(na, std::max(nb, std::max(nc, nd)));
    size_t s = std::min(na, std::min(nb, std::min(nc, nd)));
    if (l - s >= 4) { return D; }
  }

  // find last_cccc
  size_t last_cccc = 7;
  {
    std::bitset<7> h_or = ha | hb | hc | hd;
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
    size_t pa = (ha & mask).count();
    size_t pb = (hb & mask).count();
    size_t pc = (hc & mask).count();
    size_t pd = (hd & mask).count();
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
    std::bitset<7> h_and = ha & hb & hc & hd;
    if (h_and == std::bitset<7>(0b1111110)) {
      size_t n = 0;
      if (ha[0]) n++;
      if (hb[0]) n++;
      if (hc[0]) n++;
      if (hd[0]) n++;
      if (n == 1) { return C; }
    }
  }

  return D;
}

std::pair<uint64_t, double> RunFromC(double e, std::mt19937_64 &rnd, uint64_t id) {
  std::uniform_real_distribution<double> uni;
  double lambda = 1.0 / std::log(1.0 - e);

  StateN4M7 current;
  uint64_t num_c = 0;
  uint64_t t = 0;
  int64_t residual = -1;
  while (!current.All()) {
    std::array<Action,4> acts = {C, C, C, C};
    if (residual >= 0) {
      for (size_t i = 0; i < 4; i++) {
        acts[i] = CAPRI4_Action(current, i);
        if (residual >= 0) {
          if (residual == 0) {
            acts[i] = (acts[i] == C) ? D : C;
          }
          residual--;
        }
        else if (uni(rnd) < e) {
          acts[i] = (acts[i] == C) ? D : C;
        }
      }
    }
    else if (!current.Any()) {
      size_t jump = static_cast<int>( std::log(uni(rnd)) * lambda );
      size_t dt = jump / 4;
      num_c += dt * 4;
      t += dt;
      size_t tgt = jump % 4;
      acts[tgt] = D;
      for (size_t i = tgt + 1; i < 4; i++) {
        if (uni(rnd) < e) { acts[i] = D; }
      }
    }
    else if (current.State1()) { // (ccccccc, ...,ccccccd)
      size_t jump = static_cast<int>( std::log(uni(rnd)) * lambda );
      size_t dt = jump / 4;
      if (dt >= 8) {
        current.ha.reset(); current.hb.reset(); current.hc.reset(); current.hd.reset();
        num_c += (dt-1) * 4; // since d is used once for each player
        t += dt;
        size_t tgt = jump % 4;
        acts[tgt] = D;
        for (size_t i = tgt + 1; i < 4; i++) {
          if (uni(rnd) < e) { acts[i] = D; }
        }
      }
      else {
        residual = jump;
        continue;
      }
    }
    else {
      for (size_t i = 0; i < 4; i++) {
        acts[i] = CAPRI4_Action(current, i);
        if (uni(rnd) < e) { acts[i] = (acts[i] == C) ? D : C; }
      }
    }

    current.Update(acts[0], acts[1], acts[2], acts[3]);

    int n_c = 0;
    for (size_t i = 0; i < 4; i++) {
      if (acts[i] == C) n_c++;
    }
    num_c += n_c;
    t++;
  }
  return std::make_pair(t, static_cast<double>(num_c) / 4.0 / t);
}

std::pair<uint64_t, double> RunFromD(double e, std::mt19937_64 &rnd, uint64_t id) {
  std::uniform_real_distribution<double> uni;
  double lambda = 1.0 / std::log(1.0 - e);

  StateN4M7 current(0b1111111);
  uint64_t num_c = 0;
  uint64_t t = 0;
  while (current.Any()) {
    std::array<Action,4> acts = {D, D, D, D};
    if (current.All()) {
      size_t jump = static_cast<int>( std::log(uni(rnd)) * lambda );
      size_t dt = jump / 4;
      // num_c += 0;
      t += dt;
      size_t tgt = jump % 4;
      acts[tgt] = C;
      for (size_t i = tgt + 1; i < 4; i++) {
        if (uni(rnd) < e) { acts[i] = C; }
      }
    }
    else {
      for (size_t i = 0; i < 4; i++) {
        acts[i] = CAPRI4_Action(current, i);
        if (uni(rnd) < e) { acts[i] = (acts[i] == C) ? D : C; }
      }
    }

    current.Update(acts[0], acts[1], acts[2], acts[3]);

    int n_c = 0;
    for (size_t i = 0; i < 4; i++) {
      if (acts[i] == C) n_c++;
    }
    num_c += n_c;
    t++;
  }
  return std::make_pair(t, static_cast<double>(num_c) / 4.0 / t);
}

int main(int argc, char* argv[]) {
  if (argc != 4) {
    std::cerr << "invalid number of arguments" << std::endl;
    std::cerr << "usage: " << argv[0] << " <n_sample> <e> <seed>" << std::endl;
  }

  uint64_t n_samples = std::strtoull(argv[1], nullptr,0);
  double e = std::strtod(argv[2], nullptr);
  uint64_t seed_base = std::strtoull(argv[3], nullptr,0);

  std::vector<std::pair<uint64_t, double> > resultsC(n_samples), resultsD(n_samples);

#pragma omp parallel for schedule(dynamic,1)
  for (uint64_t n = 0; n < n_samples * 2; n++) {
    std::cerr << "step: " << n << std::endl;
    std::seed_seq seed = {seed_base, n};
    std::mt19937_64 rnd(seed);
    if (n % 2 == 0) {
      auto t_c1 = RunFromC(e, rnd, n);
      resultsC[n/2] = t_c1;
    }
    else {
      auto t_c2 = RunFromD(e, rnd, n);
      resultsD[n/2] = t_c2;
    }
  }

  double coop_level_c = 0.0, coop_level_d = 0.0;
  uint64_t Tc = 0, Td = 0;
  for (auto x: resultsC) {
    Tc += x.first;
    coop_level_c += x.first * x.second;
  }
  for (auto x: resultsD) {
    Td += x.first;
    coop_level_d += x.first * x.second;
  }

  std::cout << Tc << ' ' << coop_level_c / Tc << std::endl;
  std::cout << Td << ' ' << coop_level_d / Td << std::endl;
  std::cout << Tc+Td << ' ' << (coop_level_c + coop_level_d) / (Tc + Td) << std::endl;
  return 0;
}
