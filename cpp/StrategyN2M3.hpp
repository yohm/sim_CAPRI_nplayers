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
#include <bitset>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "Action.hpp"
#include "UnionFind.hpp"
#include "Partition.hpp"
#include "DirectedGraph.hpp"

#ifndef STRATEGY_N2M3_HPP
#define STRATEGY_N2M3_HPP

class StrategyN2M3;

class StateN2M3 {
 public:
  StateN2M3(Action _a_3, Action _a_2, Action _a_1, Action _b_3, Action _b_2, Action _b_1) :
      a_3(_a_3), a_2(_a_2), a_1(_a_1), b_3(_b_3), b_2(_b_2), b_1(_b_1) {};
  StateN2M3(uint64_t id) :  // upper bit [a_3,a_2,a_1,b_3,b_2,b_1] lower bit
      a_3(((id >> 5u) & 1ull) ? D : C), a_2(((id >> 4u) & 1ull) ? D : C), a_1(((id >> 3u) & 1ull) ? D : C),
      b_3(((id >> 2u) & 1ull) ? D : C), b_2(((id >> 1u) & 1ull) ? D : C), b_1(((id >> 0u) & 1ull) ? D : C) {};
  StateN2M3(const char str[6]) :
      a_3(C2A(str[0])), a_2(C2A(str[1])), a_1(C2A(str[2])),
      b_3(C2A(str[3])), b_2(C2A(str[4])), b_1(C2A(str[5])) {};
  Action a_3, a_2, a_1, b_3, b_2, b_1;

  bool operator==(const StateN2M3 &rhs) const {
    return (a_3 == rhs.a_3 && a_2 == rhs.a_2 && a_1 == rhs.a_1 && b_3 == rhs.b_3 && b_2 == rhs.b_2 && b_1 == rhs.b_1);
  }
  friend std::ostream &operator<<(std::ostream &os, const StateN2M3 &s) {
    os << s.a_3 << s.a_2 << s.a_1 << s.b_3 << s.b_2 << s.b_1;
    return os;
  };

  StateN2M3 NextState(Action act_a, Action act_b) const {
    return StateN2M3(a_2, a_1, act_a, b_2, b_1, act_b);
  };

  std::array<StateN2M3, 4> PossiblePrevStates() const {
    std::array<StateN2M3, 4> ans = {
        StateN2M3(C, a_3, a_2, C, b_3, b_2),
        StateN2M3(C, a_3, a_2, D, b_3, b_2),
        StateN2M3(D, a_3, a_2, C, b_3, b_2),
        StateN2M3(D, a_3, a_2, D, b_3, b_2)
    };
    return ans;
  }

  int RelativePayoff() const {
    if (a_1 == C && b_1 == D) { return -1; }
    else if (a_1 == D && b_1 == C) { return 1; }
    else if (a_1 == b_1) { return 0; }
    else {
      assert(false);
      return -10000;
    }
  }

  StateN2M3 SwapAB() const { return StateN2M3(b_3, b_2, b_1, a_3, a_2, a_1); } // state from B's viewpoint

  std::array<StateN2M3, 2> NoisedStates() const {
    Action a_1_n = (a_1 == C) ? D : C;
    Action b_1_n = (b_1 == C) ? D : C;
    std::array<StateN2M3, 2> ans = {StateN2M3(a_3, a_2, a_1_n, b_3, b_2, b_1), StateN2M3(a_3, a_2, a_1, b_3, b_2, b_1_n)};
    return ans;
  }
  int NumDiffInT1(const StateN2M3 &other) const {
    if (a_3 != other.a_3 || a_2 != other.a_2 || b_3 != other.b_3 || b_2 != other.b_2) {
      return -1;
    } else {
      int cnt = 0;
      if (a_1 != other.a_1) cnt++;
      if (b_1 != other.b_1) cnt++;
      return cnt;
    }
  }

  uint64_t ID() const {  // ID must be 0~63 integer. AllC: 0, AllD: 63
    uint64_t id = 0;
    if (a_3 == D) { id += 1ull << 5u; }
    if (a_2 == D) { id += 1ull << 4u; }
    if (a_1 == D) { id += 1ull << 3u; }
    if (b_3 == D) { id += 1ull << 2u; }
    if (b_2 == D) { id += 1ull << 1u; }
    if (b_1 == D) { id += 1ull << 0u; }
    return id;
  }
  bool operator<(const StateN2M3 &rhs) const {
    return (ID() < rhs.ID());
  }
};


class StrategyN2M3 {
 public:
  explicit StrategyN2M3(const std::array<Action, 64> &acts); // construct a strategy from a list of actions
  explicit StrategyN2M3(uint64_t acts);
  explicit StrategyN2M3(const char acts[64]);
  std::array<Action, 64> actions;

  std::string ToString() const;
  friend std::ostream &operator<<(std::ostream &os, const StrategyN2M3 &strategy);
  bool operator==(const StrategyN2M3 &rhs) const {
    for (size_t i = 0; i < 64; i++) { if (actions[i] != rhs.actions[i]) return false; }
    return true;
  }

  Action ActionAt(const StateN2M3 &s) const { return actions[s.ID()]; }
  void SetAction(const StateN2M3 &s, Action a) { actions[s.ID()] = a; }
  bool IsDefensible() const;  // check defensibility.
  bool IsDefensibleDFA() const; // check defensibility using DFA minimization
  // get stationary state. When coplayer is nullptr, it is set to self
  std::array<double, 64> StationaryState(double e = 0.0001, const StrategyN2M3 *coplayer = nullptr) const;
  std::array<double, 64> StationaryState2(double e = 0.0001, const StrategyN2M3 *coplayer = nullptr) const;
  // check efficiency. all actions must be fixed
  bool IsEfficient(double e = 0.00001, double th = 0.95) const { return (StationaryState(e)[0] > th); }
  bool IsEfficientTopo() const; // check efficiency using ITG
  bool IsDistinguishable(double e = 0.00001, double th = 0.95) const {
    const StrategyN2M3 allc("cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc");
    auto s = StationaryState(e, &allc);
    return (s[0] < th);
  };  // check distinguishability against AllC
  bool IsDistinguishableTopo() const; // check distinguishability using the transition graph
  DirectedGraph ITG() const;  // construct g(S,S).
  std::array<int, 64> DestsOfITG() const; // Trace g(S,S) from node i. Destination is stored in i'th element.
  int NextITGState(const StateN2M3 &s) const; // Trace the intra-transition graph by one step
  UnionFind MinimizeDFA(bool noisy = false) const;
  Partition MinimizeDFAHopcroft(bool noisy) const;
  static StrategyN2M3 CAPRI2();
  static StrategyN2M3 sCAPRI2();
  static StrategyN2M3 CAPRI() { return StrategyN2M3("cdddcdddcdcddddddcddcdddddddddddcdcdcdcdddddddddddddcdccddddddcd"); }
  static StrategyN2M3 TFT_ATFT() { return StrategyN2M3("cdcdcdcddccddccdcdcccdccdccddccdcdcdcdcddccddccdcdcccdccdccddccd"); }
  static StrategyN2M3 AON(size_t n);
 private:
  typedef std::array<std::array<int8_t, 64>, 64> d_matrix_t;
  std::vector<StateN2M3> NextPossibleStates(StateN2M3 current) const;
  bool _Equivalent(size_t i, size_t j, UnionFind &uf_0, bool noisy) const;
  typedef std::pair<size_t, int> splitter_t;
  std::array<std::set<size_t>,2> _SplitBySplitter(const Partition &partition, size_t org, const std::set<size_t> &Q, int b, bool noisy) const;
  static Action CAPRI2_Action_at(size_t i, bool is_scapri);
};

#endif //STRATEGY_N2M3_HPP

