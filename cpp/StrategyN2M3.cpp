#include <iostream>
#include <set>
#include <map>
#include "StrategyN2M3.hpp"

StrategyN2M3::StrategyN2M3(const std::array<Action, 64> &acts) : actions(acts) {}

StrategyN2M3::StrategyN2M3(const char acts[64]) {
  for (size_t i = 0; i < 64; i++) {
    actions[i] = C2A(acts[i]);
  }
}

std::vector<StateN2M3> StrategyN2M3::NextPossibleStates(StateN2M3 current) const {
  std::vector<StateN2M3> next_states;
  Action act_a = ActionAt(current);
  next_states.push_back(current.NextState(act_a, C));
  next_states.push_back(current.NextState(act_a, D));
  return std::move(next_states);
}

std::ostream &operator<<(std::ostream &os, const StrategyN2M3 &strategy) {
  for (size_t i = 0; i < 64; i++) {
    os << strategy.actions[i] << '|' << StateN2M3(i) << "  ";
    if (i % 8 == 7) { os << std::endl; }
  }
  return os;
}

std::string StrategyN2M3::ToString() const {
  char c[65];
  for (size_t i = 0; i < 64; i++) {
    c[i] = A2C(actions[i]);
  }
  c[64] = '\0';
  return std::string(c);
}

inline int8_t MIN(int8_t a, int8_t b) { return (a < b) ? a : b; }

bool StrategyN2M3::IsDefensible() const {
  const size_t N = 64;

  d_matrix_t d;

  // construct adjacency matrix
  const int INF = 32; // 32 is large enough since the path length is between -16 to 16.
  for (size_t i = 0; i < N; i++) {
    for (size_t j = 0; j < N; j++) {
      d[i][j] = INF;
    }
  }

  for (size_t i = 0; i < N; i++) {
    StateN2M3 si(i);
    std::vector<StateN2M3> sjs = NextPossibleStates(si);
    for (auto sj: sjs) {
      size_t j = sj.ID();
      d[i][j] = si.RelativePayoff();
    }
    if (d[i][i] < 0) { return false; }
  }

  for (size_t k = 0; k < N; k++) {
    for (size_t i = 0; i < N; i++) {
      for (size_t j = 0; j < N; j++) {
        d[i][j] = MIN(d[i][j], d[i][k] + d[k][j]);
      }
      if (d[i][i] < 0) { return false; }
    }
  }
  return true;
}

bool StrategyN2M3::IsDefensibleDFA() const {
  const auto autom = MinimizeDFA(false).to_map();
  std::vector<std::set<size_t> > groups;
  groups.reserve(autom.size());
  for (const auto &kv : autom) { groups.emplace_back(kv.second); }
  const size_t AN = groups.size(); // automaton size

  // initialize d_matrix
  const int INF = 32; // 32 is large enough since the path length is between -16 to 16.
  std::vector<std::vector<int> > d(AN);
  for (size_t i = 0; i < AN; i++) { d[i].resize(AN, INF); }

  auto group_index_of_state = [&groups](const StateN2M3 &ns) -> size_t {
    auto found = std::find_if(groups.cbegin(), groups.cend(), [&ns](std::set<size_t> g) {
      return g.find(ns.ID()) != g.cend();
    });
    assert(found != groups.cend());
    return std::distance(groups.cbegin(), found);
  };

  // set distance matrix
  for (size_t i = 0; i < AN; i++) {
    StateN2M3 sa = StateN2M3(*groups[i].cbegin());  // get a first state
    Action act_a = ActionAt(sa);  // A's action is same for all states in this group
    for (const auto act_b: std::array<Action, 2>({C, D})) {
      StateN2M3 ns = sa.NextState(act_a, act_b);
      size_t j = group_index_of_state(ns);
      d[i][j] = MIN(d[i][j], ns.RelativePayoff());
    }
    if (d[i][i] < 0) { return false; }
  }

  for (size_t k = 0; k < AN; k++) {
    for (size_t i = 0; i < AN; i++) {
      for (size_t j = 0; j < AN; j++) {
        d[i][j] = MIN(d[i][j], d[i][k] + d[k][j]);
      }
      if (d[i][i] < 0) { return false; }
    }
  }
  return true;
}

std::array<int, 64> StrategyN2M3::DestsOfITG() const {
  std::array<int, 64> dests = {};
  std::array<bool, 64> fixed = {false};

  for (int i = 0; i < 64; i++) {
    std::array<bool, 64> visited = {false}; // initialize by false
    visited[i] = true;
    StateN2M3 init(i);
    int next = NextITGState(init);
    while (next >= 0) {
      if (visited[next] || fixed[next]) { break; }
      visited[next] = true;
      next = NextITGState(StateN2M3(next));
    }
    int d = next;
    if (next >= 0) {
      d = fixed[next] ? dests[next] : next;
    }
    for (uint64_t j = 0; j < 64; j++) {
      if (visited[j]) {
        dests[j] = d;
        fixed[j] = true;
      }
    }
  }
  return dests;
}

int StrategyN2M3::NextITGState(const StateN2M3 &s) const {
  Action move_a = ActionAt(s);
  Action move_b = ActionAt(s.SwapAB());
  if ((move_a == C || move_a == D) && (move_b == C || move_b == D)) {
    return s.NextState(move_a, move_b).ID();
  }
  return -1;
}

std::array<double, 64> StrategyN2M3::StationaryState2(double e, const StrategyN2M3 *coplayer) const {
  if (coplayer == NULL) { coplayer = this; }
  Eigen::Matrix<double, 65, 64> A;

  for (int i = 0; i < 64; i++) {
    const StateN2M3 si(i);
    for (int j = 0; j < 64; j++) {
      // calculate transition probability from j to i
      const StateN2M3 sj(j);
      // StateN2M3 next = NextITGState(sj);
      Action act_a = ActionAt(sj);
      Action act_b = coplayer->ActionAt(sj.SwapAB());
      StateN2M3 next = sj.NextState(act_a, act_b);
      int d = next.NumDiffInT1(si);
      if (d < 0) {
        A(i, j) = 0.0;
      } else if (d == 0) {
        A(i, j) = (1.0 - e) * (1.0 - e);
      } else if (d == 1) {
        A(i, j) = (1.0 - e) * e;
      } else if (d == 2) {
        A(i, j) = e * e;
      } else {
        assert(false);
      }
    }
    A(i, i) = A(i, i) - 1.0;  // subtract unit matrix
  }
  for (int i = 0; i < 64; i++) { A(64, i) = 1.0; }  // normalization condition

  Eigen::VectorXd b(65);
  for (int i = 0; i < 64; i++) { b(i) = 0.0; }
  b(64) = 1.0;

  Eigen::VectorXd x = A.colPivHouseholderQr().solve(b);

  std::array<double, 64> ans = {0};
  for (int i = 0; i < 64; i++) {
    ans[i] = x(i);
  }
  return ans;
}

std::array<double, 64> StrategyN2M3::StationaryState(double e, const StrategyN2M3 *coplayer) const {
  if (coplayer == NULL) { coplayer = this; }

  typedef Eigen::Triplet<double> T;
  std::vector<T> tripletVec;

  for (int i = 0; i < 64; i++) {
    const StateN2M3 si(i);
    for (int j = 0; j < 64; j++) {
      // calculate transition probability from j to i
      const StateN2M3 sj(j);
      Action act_a = ActionAt(sj);
      Action act_b = coplayer->ActionAt(sj.SwapAB());
      StateN2M3 next = sj.NextState(act_a, act_b);
      int d = next.NumDiffInT1(si);
      double aij = 0.0;
      if (d < 0) {
        // A(i, j) = 0.0;
      } else if (d == 0) {
        // A(i, j) = (1.0 - e) * (1.0 - e);
        aij = (1.0 - e) * (1.0 - e);
      } else if (d == 1) {
        // A(i, j) = (1.0 - e) * e;
        aij = (1.0 - e) * e;
      } else if (d == 2) {
        // A(i, j) = e * e;
        aij = e * e;
      } else {
        assert(false);
      }
      if (i == j) { aij -= 1.0; }  // subtract unix matrix
      if (i == 63) { aij += 1.0; } // normalization condition
      if (aij != 0.0) {
        tripletVec.emplace_back(i, j, aij);
      }
    }
  }
  Eigen::SparseMatrix<double> A(64, 64);
  A.setFromTriplets(tripletVec.cbegin(), tripletVec.cend());

  Eigen::VectorXd b(64);
  for (int i = 0; i < 63; i++) { b(i) = 0.0; }
  b(63) = 1.0;

  // Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > solver;
  Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double> > solver;
  solver.compute(A);
  Eigen::VectorXd x = solver.solve(b);

  std::cerr << "#iterations:     " << solver.iterations() << std::endl;
  std::cerr << "estimated error: " << solver.error() << std::endl;

  std::array<double, 64> ans = {0};
  for (int i = 0; i < 64; i++) {
    ans[i] = x(i);
  }
  return ans;
}

DirectedGraph StrategyN2M3::ITG() const {
  DirectedGraph g(64);
  for (int i = 0; i < 64; i++) {
    StateN2M3 sa(i);
    Action act_a = ActionAt(sa);
    Action act_b = ActionAt(sa.SwapAB());
    int n = sa.NextState(act_a, act_b).ID();
    g.AddLink(i, n);
  }
  return std::move(g);
}

bool StrategyN2M3::IsEfficientTopo() const {
  if (actions[0] != C) { return false; }

  auto UpdateGn = [](DirectedGraph &gn) {
    components_t sinks = gn.SinkSCCs();
    for (const comp_t &sink: sinks) {
      for (long from: sink) {
        StateN2M3 s_from(from);
        for (const auto &s_to : s_from.NoisedStates()) {
          uint64_t to = s_to.ID();
          if (!gn.HasLink(from, to)) {
            gn.AddLink(from, to);
          }
        }
      }
    }
  };

  std::vector<int> checked(64, 0);
  checked[0] = 1;
  auto complete = [&checked]() {
    for (int i: checked) {
      if (i == 0) { return false; }
    }
    return true;
  };

  DirectedGraph gn = ITG();
  for (int n = 0; !complete(); n++) {
    if (n > 0) {
      UpdateGn(gn);
    }
    for (int i = 1; i < 64; i++) {
      if (checked[i] == 1) { continue; }
      if (gn.Reachable(i, 0)) {
        if (gn.Reachable(0, i)) {
          return false;   // inefficient
        } else {
          checked[i] = 1;
        }
      }
    }
  }

  return true;
}
bool StrategyN2M3::IsDistinguishableTopo() const {
  const StrategyN2M3 allc("cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc");
  if (actions[0] != C) { return true; }

  auto UpdateGn = [](DirectedGraph &gn) {
    components_t sinks = gn.SinkSCCs();
    for (const comp_t &sink: sinks) {
      for (long from: sink) {
        StateN2M3 s_from(from);
        for (const auto &s_to : s_from.NoisedStates()) {
          uint64_t to = s_to.ID();
          if (!gn.HasLink(from, to)) {
            gn.AddLink(from, to);
          }
        }
      }
    }
  };

  std::vector<int> checked(64, 0);
  checked[0] = 1;
  auto complete = [&checked]() {
    for (int i: checked) {
      if (i == 0) { return false; }
    }
    return true;
  };

  DirectedGraph gn(64);
  for (int i = 0; i < 64; i++) {
    StateN2M3 sa(i);
    StateN2M3 sb = sa.SwapAB();
    Action act_a = ActionAt(sa);
    Action act_b = allc.ActionAt(sb);
    assert(act_b == C);  // assert AllC
    int j = sa.NextState(act_a, act_b).ID();
    gn.AddLink(i, j);
  }

  for (int n = 0; !complete(); n++) {
    if (n > 0) {
      UpdateGn(gn);
    }
    for (int i = 1; i < 64; i++) {
      if (checked[i] == 1) { continue; }
      if (gn.Reachable(i, 0)) {
        if (gn.Reachable(0, i)) {
          return true;   // inefficient
        } else {
          checked[i] = 1;
        }
      }
    }
  }

  return false;
}

UnionFind StrategyN2M3::MinimizeDFA(bool noisy) const {
  UnionFind uf_0(64);
  // initialize grouping by the action c/d
  size_t c_rep = -1, d_rep = -1;
  for (size_t i = 0; i < 64; i++) {
    if (actions[i] == C) {
      c_rep = i;
      break;
    }
  }
  for (size_t i = 0; i < 64; i++) {
    if (actions[i] == D) {
      d_rep = i;
      break;
    }
  }
  for (size_t i = 0; i < 64; i++) {
    size_t target = (actions[i] == C) ? c_rep : d_rep;
    uf_0.merge(i, target);
  }

  while (true) {
    UnionFind uf(64);
    const auto uf_0_map = uf_0.to_map();
    for (const auto &kv : uf_0_map) { // refining a set in uf_0
      const auto &group = kv.second;
      // iterate over combinations in group
      for (auto it_i = group.cbegin(); it_i != group.cend(); it_i++) {
        auto it_j = it_i;
        it_j++;
        for (; it_j != group.cend(); it_j++) {
          if (_Equivalent(*it_i, *it_j, uf_0, noisy)) {
            uf.merge(*it_i, *it_j);
          }
        }
      }
    }
    if (uf_0_map == uf.to_map()) break;
    uf_0 = uf;
  }
  return uf_0;
}

bool StrategyN2M3::_Equivalent(size_t i, size_t j, UnionFind &uf_0, bool noisy) const {
  assert(actions[i] == actions[j]);
  Action act_a = actions[i];
  Action err_a = (act_a == C) ? D : C;
  std::array<Action, 2> acts_b = {C, D};
  for (const Action &act_b : acts_b) {
    size_t ni = StateN2M3(i).NextState(act_a, act_b).ID();
    size_t nj = StateN2M3(j).NextState(act_a, act_b).ID();
    if (uf_0.root(ni) != uf_0.root(nj)) { return false; }
    if (noisy) {
      size_t ni2 = StateN2M3(i).NextState(err_a, act_b).ID();
      size_t nj2 = StateN2M3(j).NextState(err_a, act_b).ID();
      if (uf_0.root(ni2) != uf_0.root(nj2)) { return false; }
    }
  }
  return true;
}