#include <iostream>
#include <set>
#include <map>
#include "StrategyN3M5.hpp"


const size_t StrategyN3M5::N;

StrategyN3M5::StrategyN3M5(const std::bitset<N> &acts) : actions(acts) {}

std::string StrategyN3M5::ToString() const {
  std::string s = actions.to_string();
  std::reverse(s.begin(), s.end());
  return s;
}

std::ostream &operator<<(std::ostream &os, const StrategyN3M5 &strategy) {
  for (size_t i = 0; i < StrategyN3M5::N; i++) {
    os << strategy.actions[i] << '|' << StateN3M5(i) << "  ";
    if (i % 8 == 7) { os << std::endl; }
  }
  return os;
}

std::vector<StateN3M5> StrategyN3M5::NextPossibleStates(StateN3M5 current) const {
  std::vector<StateN3M5> next_states;
  Action act_a = ActionAt(current);
  next_states.push_back(current.NextState(act_a, C, C));
  next_states.push_back(current.NextState(act_a, C, D));
  next_states.push_back(current.NextState(act_a, D, C));
  next_states.push_back(current.NextState(act_a, D, D));
  return std::move(next_states);
}

inline int16_t MIN(int16_t a, int16_t b) { return (a < b) ? a : b; }

bool StrategyN3M5::IsDefensible() const {
  typedef std::array<std::array<int16_t, N>, N> d_matrix_t;
  d_matrix_t d;

  // construct adjacency matrix
  const int INF = N/2; // N is large enough since the path length is between -N/4 to N/4.
  for (size_t i = 0; i < N; i++) {
    for (size_t j = 0; j < N; j++) {
      d[i][j] = INF;
    }
  }

  for (size_t i = 0; i < N; i++) {
    StateN3M5 si(i);
    std::vector<StateN3M5> sjs = NextPossibleStates(si);
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

bool StrategyN3M5::IsDefensibleDFA() const {
  const auto autom = MinimizeDFA(false).to_map();
  std::vector<std::set<size_t> > groups;
  groups.reserve(autom.size());
  for (const auto &kv : autom) { groups.emplace_back(kv.second); }
  const size_t AN = groups.size(); // automaton size

  // initialize d_matrix
  const int INF = N/2; // N is large enough since the path length is between -N/4 to N/4.
  std::vector<std::vector<int> > d(AN);
  for (size_t i = 0; i < AN; i++) { d[i].resize(AN, INF); }

  auto group_index_of_state = [&groups](const StateN3M5 &ns) -> size_t {
    auto found = std::find_if(groups.cbegin(), groups.cend(), [&ns](std::set<size_t> g) {
      return g.find(ns.ID()) != g.cend();
    });
    assert(found != groups.cend());
    return std::distance(groups.cbegin(), found);
  };

  // set distance matrix
  for (size_t i = 0; i < AN; i++) {
    StateN3M5 sa = StateN3M5(*groups[i].cbegin());  // get a first state
    Action act_a = ActionAt(sa);  // A's action is same for all states in this group
    for (const auto act_b: std::array<Action, 2>({C, D})) {
      for (const auto act_c: std::array<Action, 2>({C, D})) {
        StateN3M5 ns = sa.NextState(act_a, act_b, act_c);
        size_t j = group_index_of_state(ns);
        d[i][j] = MIN(d[i][j], ns.RelativePayoff());
      }
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

std::array<uint64_t , StrategyN3M5::N> StrategyN3M5::DestsOfITG() const {
  std::array<u_int64_t , N> dests = {};
  std::array<bool, N> fixed = {false};

  for (int i = 0; i < N; i++) {
    std::array<bool, N> visited = {false}; // initialize by false
    visited[i] = true;
    StateN3M5 init(i);
    int next = NextITGState(init);
    while (next >= 0) {
      if (visited[next] || fixed[next]) { break; }
      visited[next] = true;
      next = NextITGState(StateN3M5(next));
    }
    int d = fixed[next] ? dests[next] : next;
    for (uint64_t j = 0; j < N; j++) {
      if (visited[j]) {
        dests[j] = d;
        fixed[j] = true;
      }
    }
  }
  return dests;
}

uint64_t StrategyN3M5::NextITGState(const StateN3M5 &s) const {
  Action move_a = ActionAt(s);
  Action move_b = ActionAt(s.StateFromB());
  Action move_c = ActionAt(s.StateFromC());
  return s.NextState(move_a, move_b, move_c).ID();
}

/*
std::array<double, StrategyN3M5::N> StrategyN3M5::StationaryState2(double e, const StrategyN3M5 *B, const StrategyN3M5 *C) const {
  if (B == nullptr) { B = this; }
  if (C == nullptr) { C = this; }
  Eigen::Matrix<double, N, N> A;

  for (int i = 0; i < N; i++) {
    const StateN3M5 si(i);
    for (int j = 0; j < N; j++) {
      // calculate transition probability from j to i
      const StateN3M5 sj(j);
      Action act_a = ActionAt(sj);
      Action act_b = B->ActionAt(sj.StateFromB());
      Action act_c = C->ActionAt(sj.StateFromC());
      StateN3M5 next = sj.NextState(act_a, act_b, act_c);
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
  for (int i = 0; i < N; i++) { A(N-1, i) += 1.0; }  // normalization condition

  Eigen::VectorXd b(N);
  for (int i = 0; i < N-1; i++) { b(i) = 0.0; }
  b(N-1) = 1.0;

  Eigen::VectorXd x = A.colPivHouseholderQr().solve(b);

  std::array<double, N> ans = {0.0};
  for (int i = 0; i < N; i++) { ans[i] = x(i); }
  return ans;
}
 */

std::array<double, StrategyN3M5::N> StrategyN3M5::StationaryState(double e, const StrategyN3M5 *B, const StrategyN3M5 *C) const {
  if (B == nullptr) { B = this; }
  if (C == nullptr) { C = this; }

  typedef Eigen::Triplet<double> T;
  std::vector<T> tripletVec;

  for (int i = 0; i < N; i++) {
    std::cerr << "i : " << i << std::endl;
    const StateN3M5 si(i);
    for (const StateN3M5 & sj : si.PossiblePrevStates()) {
      // calculate transition probability from j to i
      size_t j = sj.ID();
      Action act_a = ActionAt(sj);
      Action act_b = B->ActionAt(sj.StateFromB());
      Action act_c = C->ActionAt(sj.StateFromC());
      StateN3M5 next = sj.NextState(act_a, act_b, act_c);
      int d = next.NumDiffInT1(si);
      double aij = 0.0;
      if (d < 0) {
        // A(i, j) = 0.0;
      } else if (d == 0) {
        // A(i, j) = (1.0 - e) * (1.0 - e);
        aij = (1.0 - e) * (1.0 - e) * (1.0 - e);
      } else if (d == 1) {
        // A(i, j) = (1.0 - e) * e;
        aij = (1.0 - e) * (1.0 - e) * e;
      } else if (d == 2) {
        // A(i, j) = e * e;
        aij = (1.0 - e) * e * e;
      } else if (d == 3) {
        aij = e * e * e;
      } else {
        assert(false);
      }
      if (aij != 0.0) {
        tripletVec.emplace_back(i, j, aij);
      }
    }
  }
  Eigen::SparseMatrix<double> A(N, N);
  A.setFromTriplets(tripletVec.cbegin(), tripletVec.cend());

  // subtract unit matrix & normalization condition
  std::vector<T> iVec;
  for (int i = 0; i < N-1; i++) { iVec.emplace_back(i, i, -1.0); }
  for (int i = 0; i < N-1; i++) { iVec.emplace_back(N-1, i, 1.0); }
  Eigen::SparseMatrix<double> I(N, N);
  I.setFromTriplets(iVec.cbegin(), iVec.cend());
  std::cerr << iVec.rbegin()->col() << ' ' << iVec.rbegin()->row() << ' ' << iVec.rbegin()->value() << std::endl;
  A = A + I;

  Eigen::VectorXd b = Eigen::VectorXd::Zero(N);
  b(N-1) = 1.0;
  std::cerr << b << std::endl;

  // Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > solver;
  Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double> > solver;
  solver.compute(A);
  Eigen::VectorXd x = solver.solve(b);

  std::cerr << "#iterations:     " << solver.iterations() << std::endl;
  std::cerr << "estimated error: " << solver.error() << std::endl;
  std::cerr << "x: " << x[0] << ' ' << x[32767] << std::endl;

  std::array<double, N> ans = {0};
  for (int i = 0; i < N; i++) { ans[i] = x[i]; }
  return ans;
}

DirectedGraph StrategyN3M5::ITG() const {
  DirectedGraph g(N);
  for (int i = 0; i < N; i++) {
    StateN3M5 sa(i);
    StateN3M5 sb = sa.StateFromB();
    StateN3M5 sc = sa.StateFromC();
    Action act_a = ActionAt(sa);
    Action act_b = ActionAt(sb);
    Action act_c = ActionAt(sc);
    int n = sa.NextState(act_a, act_b, act_c).ID();
    g.AddLink(i, n);
  }

  return std::move(g);
}

bool StrategyN3M5::IsEfficientTopo() const {
  if (actions[0]) { return false; }

  auto UpdateGn = [](DirectedGraph &gn) {
    components_t sinks = gn.SinkSCCs();
    for (const comp_t &sink: sinks) {
      for (uint64_t from: sink) {
        StateN3M5 s_from(from);
        for (const StateN3M5 &s_to: s_from.NoisedStates()) {
          uint64_t to = s_to.ID();
          if (!gn.HasLink(from, to)) {
            gn.AddLink(from, to);
          }
        }
      }
    }
  };

  std::vector<int> checked(N, 0);
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
    for (int i = 1; i < N; i++) {
      if(i%1000 == 0) {std::cerr << "checking efficiency (n,i) : " << n << ", " << i << std::endl;}
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
bool StrategyN3M5::IsDistinguishableTopo() const {
  const StrategyN3M5 allc(std::bitset<N>(0ull));
  if (actions[0]) { return true; }

  auto UpdateGn = [](DirectedGraph &gn) {
    components_t sinks = gn.SinkSCCs();
    for (const comp_t &sink: sinks) {
      for (u_int64_t from: sink) {
        StateN3M5 s_from(from);
        for (const StateN3M5 &s_to: s_from.NoisedStates()) {
          uint64_t to = s_to.ID();
          if (!gn.HasLink(from, to)) {
            gn.AddLink(from, to);
          }
        }
      }
    }
  };

  std::vector<int> checked(N, 0);
  checked[0] = 1;
  auto complete = [&checked]() {
    for (int i: checked) {
      if (i == 0) { return false; }
    }
    return true;
  };

  DirectedGraph gn(N);
  for (int i = 0; i < N; i++) {
    StateN3M5 sa(i);
    StateN3M5 sb = sa.StateFromB();
    StateN3M5 sc = sa.StateFromC();
    Action act_a = ActionAt(sa);
    Action act_b = allc.ActionAt(sb);
    Action act_c = allc.ActionAt(sc);
    assert(act_b == C && act_c == C);  // assert AllC
    int j = sa.NextState(act_a, act_b, act_c).ID();
    gn.AddLink(i, j);
  }

  for (int n = 0; !complete(); n++) {
    if (n > 0) {
      UpdateGn(gn);
    }
    for (int i = 1; i < N; i++) {
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

UnionFind StrategyN3M5::MinimizeDFA(bool noisy) const {
  UnionFind uf_0(N);
  // initialize grouping by the action c/d
  long c_rep = -1, d_rep = -1;
  for (size_t i = 0; i < N; i++) {
    if (actions[i] == false) {
      c_rep = i;
      break;
    }
  }
  for (size_t i = 0; i < N; i++) {
    if (actions[i] == true) {
      d_rep = i;
      break;
    }
  }
  for (size_t i = 0; i < N; i++) {
    long target = (actions[i] == false) ? c_rep : d_rep;
    uf_0.merge(i, target);
  }

  while (true) {
    UnionFind uf(N);
    const auto uf_0_map = uf_0.to_map();
    for (const auto &kv : uf_0_map) { // refining a set in uf_0
      const auto &group = kv.second;
      // iterate over combinations in group
      size_t idx = 0, n = group.size();
      for (auto it_i = group.cbegin(); it_i != group.cend(); it_i++, idx++) {
        auto it_j = it_i;
        std::cerr << "idx / n : " << idx << " / " << n << std::endl;
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

bool StrategyN3M5::_Equivalent(size_t i, size_t j, UnionFind &uf_0, bool noisy) const {
  assert(actions[i] == actions[j]);
  Action act_a = ActionAt(i);
  Action err_a = (act_a == C) ? D : C;
  std::array<Action, 2> acts_b = {C, D};
  std::array<Action, 2> acts_c = {C, D};
  for (const Action &act_b : acts_b) {
    for (const Action &act_c : acts_c) {
      size_t ni = StateN3M5(i).NextState(act_a, act_b, act_c).ID();
      size_t nj = StateN3M5(j).NextState(act_a, act_b, act_c).ID();
      if (uf_0.root(ni) != uf_0.root(nj)) { return false; }
      if (noisy) {
        size_t ni2 = StateN3M5(i).NextState(err_a, act_b, act_c).ID();
        size_t nj2 = StateN3M5(j).NextState(err_a, act_b, act_c).ID();
        if (uf_0.root(ni2) != uf_0.root(nj2)) { return false; }
      }
    }
  }
  return true;
}
