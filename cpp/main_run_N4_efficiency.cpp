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
  typedef std::bitset<7> BT;
  StateN4M7(uint64_t i = 0) : ha(i), hb(i), hc(i), hd(i) {};
  // StateN4M7(BT _ha, BT _hb, BT _hc, BT _hd) : ha(_ha), hb(_hb), hc(_hc), hd(_hd) {};
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
      if (n == 1) {std::cerr << "R" << std::endl; return C; }
    }
  }

  return D;
}

void Run(double e, std::mt19937_64 &rnd, uint64_t id) {
  std::uniform_real_distribution<double> uni;
  double lambda = 1.0 / std::log(1.0 - e);

  StateN4M7 current;
  int64_t num_c = 0;
  int64_t t = 0;
  while (!current.All()) {
    std::array<Action,4> acts = {C, C, C, C};
    if (!current.Any()) {
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
  std::cout << id << ' ' << t << ' ' << static_cast<double>(num_c) / 4.0 / (t+1) << std::endl;
}

int main(int argc, char* argv[]) {
  if (argc != 4) {
    std::cerr << "invalid number of arguments" << std::endl;
    std::cerr << "usage: " << argv[0] << " <n_sample> <e> <seed>" << std::endl;
  }

  uint64_t n_samples = std::strtoull(argv[1], nullptr,0);
  double e = std::strtod(argv[2], nullptr);
  uint64_t seed_base = std::strtoull(argv[3], nullptr,0);

  for (uint64_t n = 0; n < n_samples; n++) {
    std::seed_seq seed = {seed_base, n};
    std::mt19937_64 rnd(seed);
    Run(e, rnd, n);

  }
  return 0;
}
