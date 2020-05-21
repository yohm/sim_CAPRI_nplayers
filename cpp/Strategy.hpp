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
#include "Action.hpp"
#include "DirectedGraph.hpp"

#ifndef STRATEGY_HPP
#define STRATEGY_HPP

class Strategy;

class State {
 public:
  State(Action _a_3, Action _a_2, Action _a_1, Action _b_3, Action _b_2, Action _b_1) :
      a_3(_a_3), a_2(_a_2), a_1(_a_1), b_3(_b_3), b_2(_b_2), b_1(_b_1) { assert(AllCorD()); };
  State(uint64_t id) :  // upper bit [a_3,a_2,a_1,b_3,b_2,b_1] lower bit
      a_3(((id >> 5) & 1) ? D : C), a_2(((id >> 4) & 1) ? D : C), a_1(((id >> 3) & 1) ? D : C),
      b_3(((id >> 2) & 1) ? D : C), b_2(((id >> 1) & 1) ? D : C), b_1(((id >> 0) & 1) ? D : C) { assert(AllCorD()); };
  State(const char str[6]) :
      a_3(C2A(str[0])), a_2(C2A(str[1])), a_1(C2A(str[2])),
      b_3(C2A(str[3])), b_2(C2A(str[4])), b_1(C2A(str[5])) { assert(AllCorD()); };
  const Action a_3, a_2, a_1, b_3, b_2, b_1;
  bool AllCorD() const {
    return (
        (a_3 == D || a_3 == C) && (a_2 == D || a_2 == C) && (a_1 == D || a_1 == C) &&
            (b_3 == D || b_3 == C) && (b_2 == D || b_2 == C) && (b_1 == D || b_1 == C)
    );
  }

  bool operator==(const State &rhs) const {
    return (a_3 == rhs.a_3 && a_2 == rhs.a_2 && a_1 == rhs.a_1 && b_3 == rhs.b_3 && b_2 == rhs.b_2 && b_1 == rhs.b_1);
  }
  friend std::ostream &operator<<(std::ostream &os, const State &s) {
    os << s.a_3 << s.a_2 << s.a_1 << s.b_3 << s.b_2 << s.b_1;
    return os;
  };

  State NextState(Action act_a, Action act_b) const {
    assert(act_a == C || act_a == D);
    assert(act_b == C || act_b == D);
    return State(a_2, a_1, act_a, b_2, b_1, act_b);
  };

  std::array<State, 4> PossiblePrevStates() const {
    std::array<State, 4> ans = {
        State(C, a_3, a_2, C, b_3, b_2),
        State(C, a_3, a_2, D, b_3, b_2),
        State(D, a_3, a_2, C, b_3, b_2),
        State(D, a_3, a_2, D, b_3, b_2)
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

  State SwapAB() const { return State(b_3, b_2, b_1, a_3, a_2, a_1); } // state from B's viewpoint

  std::array<State, 2> NoisedStates() const {
    Action a_1_n = (a_1 == C) ? D : C;
    Action b_1_n = (b_1 == C) ? D : C;
    std::array<State, 2> ans = {State(a_3, a_2, a_1_n, b_3, b_2, b_1), State(a_3, a_2, a_1, b_3, b_2, b_1_n)};
    return ans;
  }
  int NumDiffInT1(const State &other) const {
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
    if (a_3 == D) { id += 1 << 5; }
    if (a_2 == D) { id += 1 << 4; }
    if (a_1 == D) { id += 1 << 3; }
    if (b_3 == D) { id += 1 << 2; }
    if (b_2 == D) { id += 1 << 1; }
    if (b_1 == D) { id += 1 << 0; }
    return id;
  }
  bool operator<(const State &rhs) const {
    return (ID() < rhs.ID());
  }

};

class UnionFind {
 public:
  UnionFind(size_t n) : parent(n) {
    for (size_t i = 0; i < n; i++) { parent[i] = i; }
  }
  size_t root(size_t i) {
    if (parent[i] != i) {
      size_t r = root(parent[i]);
      parent[i] = r;
    }
    return parent[i];
  }
  bool merge(size_t i, size_t j) {
    size_t ri = root(i);
    size_t rj = root(j);
    if (ri == rj) return false;
    else if (ri > rj) { parent[ri] = rj; }
    else if (ri < rj) { parent[rj] = ri; }
    return true;
  }
  std::map<size_t, std::set<size_t> > to_map() {
    std::map<size_t, std::set<size_t> > m;
    for (size_t i = 0; i < parent.size(); i++) {
      size_t r = root(i);
      m[r].insert(i);
    }
    return std::move(m);
  }
 private:
  std::vector<size_t> parent;
};

class Strategy {
 public:
  Strategy(const std::array<Action, 64> &acts); // construct a strategy from a list of actions
  Strategy(const char acts[64]);
  std::array<Action, 64> actions;

  std::string ToString() const;
  friend std::ostream &operator<<(std::ostream &os, const Strategy &strategy);
  bool operator==(const Strategy &rhs) const {
    for (size_t i = 0; i < 64; i++) { if (actions[i] != rhs.actions[i]) return false; }
    return true;
  }

  Action ActionAt(const State &s) const { return actions[s.ID()]; }
  void SetAction(const State &s, Action a) { actions[s.ID()] = a; }
  bool IsDefensible() const;  // check defensibility.
  // get stationary state. When coplayer is nullptr, it is set to self
  std::array<double, 64> StationaryState(double e = 0.0001, const Strategy *coplayer = nullptr) const;
  // check efficiency. all actions must be fixed
  bool IsEfficient(double e = 0.00001, double th = 0.95) const { return (StationaryState(e)[0] > th); }
  bool IsEfficientTopo() const; // check efficiency using ITG
  bool IsDistinguishable(double e = 0.00001, double th = 0.95) const {
    const Strategy allc("cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc");
    return (StationaryState(e, &allc)[0] < th);
  };  // check distinguishability against AllC
  bool IsDistinguishableTopo() const; // check distinguishability using the transition graph
  DirectedGraph ITG() const;  // construct g(S,S).
  std::array<int, 64> DestsOfITG() const; // Trace g(S,S) from node i. Destination is stored in i'th element.
  int NextITGState(const State &s) const; // Trace the intra-transition graph by one step
  UnionFind MinimizeDFA(bool noisy = false) const;
 private:
  typedef std::array<std::array<int8_t, 64>, 64> d_matrix_t;
  std::vector<State> NextPossibleStates(State current) const;
  bool _Equivalent(size_t i, size_t j, UnionFind& uf_0, bool noisy) const;
};

#endif //STRATEGY_HPP

