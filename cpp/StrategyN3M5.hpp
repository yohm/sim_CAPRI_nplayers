#include <string>
#include <array>
#include <vector>
#include <set>
#include <map>
#include <queue>
#include <stack>
#include <sstream>
#include <cstdint>
#include <ostream>
#include <algorithm>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "Action.hpp"
#include "UnionFind.hpp"
#include "DirectedGraph.hpp"

#ifndef STRATEGY_N3M5_HPP
#define STRATEGY_N3M5_HPP

class StrategyN3M5;

typedef std::bitset<5> HistoM5;

class StateN3M5 {
 public:
  StateN3M5(const HistoM5 &a_histo, const HistoM5 &b_histo, const HistoM5 &c_histo) :
      ha(a_histo), hb(b_histo), hc(c_histo) {};
  StateN3M5(const std::string &str) { // format "ccccc_dddcc_cdccd" (a0,a1,a2,...c3,c4) last action comes first
    assert(str.size() == 17 && str[5] == '_' && str[11] == '_');
    std::string sa = str.substr(0, 5), sb = str.substr(6, 5), sc = str.substr(12, 5);
    auto s_to_binary = [](std::string s)->HistoM5 {
      std::reverse(s.begin(), s.end());
      std::replace(s.begin(), s.end(), 'c', '0');
      std::replace(s.begin(), s.end(), 'd', '1');
      return HistoM5(s);
    };
    ha = s_to_binary(sa); hb = s_to_binary(sb); hc = s_to_binary(sc);
  }
  StateN3M5(uint64_t i) { // (c4,c3,c2,....,a1,a0) the lowest bit is a0
    uint64_t mask = 31ull;
    ha = HistoM5(i & mask); hb = HistoM5((i>>5ul)&mask); hc = HistoM5((i>>10)&mask);
  }
  HistoM5 ha, hb, hc;

  bool operator==(const StateN3M5 &rhs) const {
    return (ha == rhs.ha && hb == rhs.hb && hc == rhs.hc);
  }
  friend std::ostream &operator<<(std::ostream &os, const StateN3M5 &s) {
    auto to_s = [](const HistoM5 &h)->std::string {
      std::string s = h.to_string();
      std::reverse(s.begin(), s.end());
      std::replace(s.begin(), s.end(), '0', 'c');
      std::replace(s.begin(), s.end(), '1', 'd');
      return s;
    };
    os << to_s(s.ha) << '_' << to_s(s.hb) << '_' << to_s(s.hc);
    return os;
  };

  StateN3M5 NextState(Action act_a, Action act_b, Action act_c) const {
    HistoM5 nha = ha << 1;
    HistoM5 nhb = hb << 1;
    HistoM5 nhc = hc << 1;
    if (act_a == D) { nha.set(0); }
    if (act_b == D) { nhb.set(0); }
    if (act_c == D) { nhc.set(0); }
    return StateN3M5(nha, nhb, nhc);
  };

  std::vector<StateN3M5> PossiblePrevStates() const {
    std::vector<StateN3M5> ans;
    for (size_t i = 0; i < 2; i++) {
      for (size_t j = 0; j < 2; j++) {
        for (size_t k = 0; k < 2; k++) {
          HistoM5 lha = ha >> 1;
          HistoM5 lhb = hb >> 1;
          HistoM5 lhc = hc >> 1;
          if (i == 1) { lha.set(4); }
          if (j == 1) { lhb.set(4); }
          if (k == 1) { lhc.set(4); }
          ans.emplace_back(lha, lhb, lhc);
        }
      }
    }
    return std::move(ans);
  }

  int RelativePayoff(bool against_B = true) const {
    bool a0 = ha[0];
    bool b0 = against_B ? hb[0] : hc[0];
    if (a0 == false && b0 == true) { return -1; } // C,D
    else if (a0 == true && b0 == false) { return 1; } // D,C
    else if (a0 == b0) { return 0; }  // C,C or D,D
    else {
      assert(false);
      return -10000;
    }
  }

  StateN3M5 StateFromB() const { return StateN3M5(hb, hc, ha); } // state from B's viewpoint
  StateN3M5 StateFromC() const { return StateN3M5(hc, ha, hb); } // state from C's viewpoint

  std::array<StateN3M5, 3> NoisedStates() const {
    HistoM5 ha_n = ha, hb_n = hb, hc_n = hc;
    ha_n ^= HistoM5(1ull);
    hb_n ^= HistoM5(1ull);
    hc_n ^= HistoM5(1ull);
    std::array<StateN3M5, 3> ans = {StateN3M5(ha_n, hb, hc), StateN3M5(ha, hb_n, hc), StateN3M5(ha, hb, hc_n)};
    return ans;
  }

  int NumDiffInT1(const StateN3M5 &other) const {
    const auto b1 = ToBits();
    const auto b2 = other.ToBits();

    std::bitset<15> mask("111101111011110");
    if( (b1 & mask) != (b2 & mask) ) { // inconsistent bit is found
      return -1;
    } else {
      return ((b1 & ~mask) ^ (b2 & ~mask)).count(); // number of different bits
    }
  }

  std::string ToString() const {
    std::ostringstream os;
    os << *this;
    return os.str();
  }

  std::bitset<15> ToBits() const {
    std::bitset<15> bits(0ull);
    bits ^= ha.to_ullong();
    bits ^= (hb.to_ullong() << 5);
    bits ^= (hc.to_ullong() <<10);
    return bits;
  }

  uint64_t ID() const {
    return ToBits().to_ullong();
  }

  bool operator<(const StateN3M5 &rhs) const {
    return (ID() < rhs.ID());
  }
};


class StrategyN3M5 {
 public:
  static const size_t N = 1ull << 15ull; // == 32768
  explicit StrategyN3M5(const std::bitset<N> &actions); // construct a strategy from a list of actions. 0=>c,1=>d

  std::string ToString() const;
  friend std::ostream &operator<<(std::ostream &os, const StrategyN3M5 &strategy);
  bool operator==(const StrategyN3M5 &rhs) const {
    for (size_t i = 0; i < 64; i++) { if (actions[i] != rhs.actions[i]) return false; }
    return true;
  }

  Action ActionAt(const StateN3M5 &s) const { return actions[s.ID()] ? D : C; }
  bool IsDefensible() const;  // check defensibility. Not computationally feasible.
  bool IsDefensibleDFA(); // check defensibility using DFA minimization
  // get stationary state. When coplayer is nullptr, it is set to self
  std::array<double, N> StationaryState(double e = 0.0001, const StrategyN3M5 *B = nullptr, const StrategyN3M5 *C = nullptr) const;
  // std::array<double, N> StationaryState2(double e = 0.0001, const StrategyN3M5 *B = nullptr, const StrategyN3M5 *C = nullptr) const;
  // check efficiency. all actions must be fixed
  bool IsEfficient(double e = 0.00001, double th = 0.95) const { return (StationaryState(e)[0] > th); }
  bool IsEfficientTopo() const; // check efficiency using ITG
  bool IsDistinguishable(double e = 0.00001, double th = 0.95) const {
    const StrategyN3M5 allc(std::bitset<N>(0ull));
    return (StationaryState(e, &allc, &allc)[0] < th);
  };  // check distinguishability against AllC
  bool IsDistinguishableTopo() const; // check distinguishability using the transition graph
  DirectedGraph ITG() const;  // construct g(S,S).
  std::array<uint64_t , StrategyN3M5::N> DestsOfITG() const; // Trace g(S,S) from node i. Destination is stored in i'th element.
  uint64_t NextITGState(const StateN3M5 &s) const; // Trace the intra-transition graph by one step
  UnionFind MinimizeDFA(bool noisy = false);
 private:
  const std::bitset<N> actions;
  UnionFind min_auto_cache[2];  // cache of simplified and full automaton
  std::vector<StateN3M5> NextPossibleStates(StateN3M5 current) const;
  bool _Equivalent(size_t i, size_t j, UnionFind &uf_0, bool noisy) const;
};


#endif //STRATEGY_N3M5_HPP
