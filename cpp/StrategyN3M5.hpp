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

typedef std::array<Action,5> HistoM5;

class StateN3M5 {
 public:
  StateN3M5(const HistoM5 &a_histo, const HistoM5 &b_histo, const HistoM5 &c_histo) :
      ha(a_histo), hb(b_histo), hc(c_histo) {};
  StateN3M5(const std::string &str) :  // format "ccccc_dddcc_cdccd" last action comes first
      ha({C2A(str[0]), C2A(str[1]), C2A(str[2]), C2A(str[3]), C2A(str[4])}),
      hb({C2A(str[6]), C2A(str[7]), C2A(str[8]), C2A(str[9]), C2A(str[10])}),
      hc({C2A(str[12]), C2A(str[13]), C2A(str[14]), C2A(str[15]), C2A(str[16])})
      { assert(str.size() == 17 && str[5] == '_' && str[11] == '_'); }
  const HistoM5 ha, hb, hc;

  bool operator==(const StateN3M5 &rhs) const {
    return (ha == rhs.ha && hb == rhs.hb && hc == rhs.hc);
  }
  friend std::ostream &operator<<(std::ostream &os, const StateN3M5 &s) {
    os << s.ha[0] << s.ha[1] << s.ha[2] << s.ha[3] << s.ha[4] << '_';
    os << s.hb[0] << s.hb[1] << s.hb[2] << s.hb[3] << s.hb[4] << '_';
    os << s.hc[0] << s.hc[1] << s.hc[2] << s.hc[3] << s.hc[4];
    return os;
  };

  StateN3M5 NextState(Action act_a, Action act_b, Action act_c) const {
    HistoM5 nha = {act_a, ha[0], ha[1], ha[2], ha[3]};
    HistoM5 nhb = {act_b, hb[0], hb[1], hb[2], hb[3]};
    HistoM5 nhc = {act_c, hc[0], hc[1], hc[2], hc[3]};
    return StateN3M5(nha, nhb, nhc);
  };

  std::vector<StateN3M5> PossiblePrevStates() const {
    std::vector<StateN3M5> ans;

    for (size_t i = 0; i < 2; i++) {
      Action act_a = (i == 0) ? C : D;
      for (size_t j = 0; j < 2; j++) {
        Action act_b = (j == 0) ? C : D;
        for (size_t k = 0; k < 2; k++) {
          Action act_c = (k == 0) ? C : D;

          HistoM5 lha = {ha[1], ha[2], ha[3], ha[4], act_a};
          HistoM5 lhb = {hb[1], hb[2], hb[3], hb[4], act_b};
          HistoM5 lhc = {hc[1], hc[2], hc[3], hc[4], act_c};
          ans.emplace_back(lha, lhb, lhc);
        }
      }
    }
    return std::move(ans);
  }

  int RelativePayoff(bool against_B = true) const {
    Action a0 = ha[0];
    Action b0 = against_B ? hb[0] : hc[0];
    if (a0 == C && b0 == D) { return -1; }
    else if (a0 == D && b0 == C) { return 1; }
    else if (a0 == b0) { return 0; }
    else {
      assert(false);
      return -10000;
    }
  }

  StateN3M5 StateFromB() const { return StateN3M5(hb, hc, ha); } // state from B's viewpoint
  StateN3M5 StateFromC() const { return StateN3M5(hc, ha, hb); } // state from C's viewpoint

  std::array<StateN3M5, 3> NoisedStates() const {
    HistoM5 ha_n = ha, hb_n = hb, hc_n = hc;
    ha_n[0] = (ha_n[0] == C) ? D : C;
    hb_n[0] = (hb_n[0] == C) ? D : C;
    hc_n[0] = (hc_n[0] == C) ? D : C;
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
    for (size_t i = 0; i < 5; i++) { if(ha[i]==D) bits.set(i); }
    for (size_t i = 0; i < 5; i++) { if(hb[i]==D) bits.set(i+5); }
    for (size_t i = 0; i < 5; i++) { if(hc[i]==D) bits.set(i+10); }
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
  explicit StrategyN3M5(const std::array<Action, 64> &acts); // construct a strategy from a list of actions
  explicit StrategyN3M5(const char acts[64]);
  std::array<Action, 64> actions;

  std::string ToString() const;
  friend std::ostream &operator<<(std::ostream &os, const StrategyN3M5 &strategy);
  bool operator==(const StrategyN3M5 &rhs) const {
    for (size_t i = 0; i < 64; i++) { if (actions[i] != rhs.actions[i]) return false; }
    return true;
  }

  Action ActionAt(const StateN3M5 &s) const { return actions[s.ID()]; }
  void SetAction(const StateN3M5 &s, Action a) { actions[s.ID()] = a; }
  bool IsDefensible() const;  // check defensibility.
  bool IsDefensibleDFA() const; // check defensibility using DFA minimization
  // get stationary state. When coplayer is nullptr, it is set to self
  std::array<double, 64> StationaryState(double e = 0.0001, const StrategyN3M5 *coplayer = nullptr) const;
  std::array<double, 64> StationaryState2(double e = 0.0001, const StrategyN3M5 *coplayer = nullptr) const;
  // check efficiency. all actions must be fixed
  bool IsEfficient(double e = 0.00001, double th = 0.95) const { return (StationaryState(e)[0] > th); }
  bool IsEfficientTopo() const; // check efficiency using ITG
  bool IsDistinguishable(double e = 0.00001, double th = 0.95) const {
    const StrategyN3M5 allc("cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc");
    return (StationaryState(e, &allc)[0] < th);
  };  // check distinguishability against AllC
  bool IsDistinguishableTopo() const; // check distinguishability using the transition graph
  DirectedGraph ITG() const;  // construct g(S,S).
  std::array<int, 64> DestsOfITG() const; // Trace g(S,S) from node i. Destination is stored in i'th element.
  int NextITGState(const StateN3M5 &s) const; // Trace the intra-transition graph by one step
  UnionFind MinimizeDFA(bool noisy = false) const;
 private:
  typedef std::array<std::array<int8_t, 64>, 64> d_matrix_t;
  std::vector<StateN3M5> NextPossibleStates(StateN3M5 current) const;
  bool _Equivalent(size_t i, size_t j, UnionFind &uf_0, bool noisy) const;
};

#endif //STRATEGY_N3M5_HPP

