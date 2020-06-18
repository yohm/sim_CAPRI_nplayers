#include <iostream>
#include <set>
#include <map>
#include "StrategyN3M5.hpp"


const size_t StrategyN3M5::N;

StrategyN3M5::StrategyN3M5(const std::bitset<N> &acts) : actions(acts) {
}

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
  const auto autom = MinimizeDFAHopcroft(false).to_map();
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

  // debug print
  auto print_d = [&d,AN]() {
    for (size_t i = 0; i < AN; i++) {
      for (size_t j = 0; j < AN; j++) {
        std::cerr << d[i][j] << ' ';
      }
      std::cerr << std::endl;
    }
  };
  // print_d();

  for (size_t k = 0; k < AN; k++) {
    for (size_t i = 0; i < AN; i++) {
      for (size_t j = 0; j < AN; j++) {
        d[i][j] = std::min(d[i][j], d[i][k] + d[k][j]);
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
    uint64_t next = NextITGState(init);
    while (true) {
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

std::vector<uint64_t> StrategyN3M5::TraceStates(uint64_t start, const StrategyN3M5 *B, const StrategyN3M5 *C) const {
  std::vector<uint64_t> trace;

  if (B == nullptr) { B = this; }
  if (C == nullptr) { C = this; }

  auto next_state = [this,B,C](uint64_t i) -> uint64_t {
    auto s = StateN3M5(i);
    Action a = this->ActionAt(s);
    Action b = B->ActionAt(s.StateFromB());
    Action c = C->ActionAt(s.StateFromC());
    return s.NextState(a, b, c).ID();
  };

  std::array<bool, N> visited = {false}; // initialize by false
  visited[start] = true;
  trace.push_back(start);
  uint64_t c = next_state(start);
  while (true) {
    if (visited[c]) { break; }
    visited[c] = true;
    trace.push_back(c);
    c = next_state(c);
  }
  return std::move(trace);
}

uint64_t StrategyN3M5::NextITGState(const StateN3M5 &s) const {
  Action move_a = ActionAt(s);
  Action move_b = ActionAt(s.StateFromB());
  Action move_c = ActionAt(s.StateFromC());
  return s.NextState(move_a, move_b, move_c).ID();
}

std::array<double, StrategyN3M5::N> StrategyN3M5::StationaryState(double e, const StrategyN3M5 *B, const StrategyN3M5 *C) const {
  std::cerr << "calculating stationary state" << std::endl;
  if (B == nullptr) { B = this; }
  if (C == nullptr) { C = this; }

  typedef Eigen::Triplet<double> T;
  std::vector<T> tripletVec;

  for (int i = 0; i < N; i++) {
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
  A = A + I;
  std::cerr << "  transition matrix has been created" << std::endl;

  Eigen::VectorXd b = Eigen::VectorXd::Zero(N);
  b(N-1) = 1.0;

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
      if(i%10000 == 0) {std::cerr << "checking efficiency (n,i) : " << n << ", " << i << std::endl;}
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
        if(idx % 100 == 0) { std::cerr << "in DFA minimization, idx / n : " << idx << " / " << n << std::endl; }
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

Partition StrategyN3M5::MinimizeDFAHopcroft(bool noisy) const {
  Partition partition(StrategyN3M5::N);

  // initialize partition by the action c/d
  std::set<size_t> c_set, d_set;
  for (size_t i = 0; i < StrategyN3M5::N; i++) {
    if (actions[i] == C) { c_set.insert(i); }
    else { d_set.insert(i); }
  }
  if (c_set.empty() || d_set.empty()) { return std::move(partition); }
  partition.split(0, c_set);
  size_t smaller = (c_set.size() < d_set.size()) ? (*c_set.begin()) : (*d_set.begin());

  std::set<splitter_t> waiting;  // initialize waiting set
  int input_size = 4;  // 0: cc, 1: cd, 2: dc, 3: dd
  if (noisy) { input_size = 8; } // 0: ccc, 1: ccd, 2: cdc, 3: cdd, 4: dcc, 5: dcd, 6: ddc, 7: ddd
  for(int b = 0; b < input_size; b++) {
    waiting.insert({smaller, b});
  }

  while( !waiting.empty() ) {
    splitter_t splitter = *waiting.cbegin();  // take a splitter from a waiting list
    const std::set<size_t> Q = partition.group(splitter.first);  // copy splitter because this must remain same during this iteration
    const int b = splitter.second;
    waiting.erase(waiting.begin());
    auto groups = partition.group_ids();
    for(size_t p: groups) {  // for each P in partition
      // check if group i is splittable or not
      auto p1p2 = _SplitBySplitter(partition, p, Q, b, noisy);
      const std::set<size_t> &p1 = p1p2.at(0), &p2 = p1p2.at(1);
      if (p1.empty() || p2.empty() ) { continue; }  // P is not split by this splitter
      partition.split(p, p1);
      for (int b = 0; b < input_size; b++) {
        auto found = waiting.find(splitter_t(p, b) );
        if (found != waiting.end()) {
          // replace (P, b) by (P1, b) and (P2, b) in W
          waiting.erase(found);
          waiting.insert({*p1.cbegin(), b});
          waiting.insert({*p2.cbegin(), b});
        }
        else {
          if (p1.size() < p2.size()) {
            waiting.insert({*p1.cbegin(), b});
          }
          else {
            waiting.insert({*p2.cbegin(), b});
          }
        }
      }
    }
  }
  return std::move(partition);
}

std::array<std::set<size_t>,2> StrategyN3M5::_SplitBySplitter(const Partition &partition, size_t p, const std::set<size_t> &Q, int b, bool noisy) const {
  const std::set<size_t> &P = partition.group(p);
  Action act_c = (b & 1) ? D : C;
  Action act_b = (b & 2) ? D : C;
  // get members of P which go to a member of Q by b
  std::set<size_t> P1;
  for (size_t si: P) {
    StateN3M5 s(si);
    Action act_a = ActionAt(s);
    if (noisy) { act_a = (b & 4) ? D : C; }
    size_t next = s.NextState(act_a, act_b, act_c).ID();
    if (Q.find(next) != Q.end()) {
      P1.insert(si);
    }
  }
  std::set<size_t> P2;
  std::set_difference(P.begin(), P.end(), P1.begin(), P1.end(), std::inserter(P2, P2.end()));
  return {P1, P2};
}


StrategyN3M5 StrategyN3M5::AllC() {
  std::bitset<N> allc_b(0ull);
  return std::move(StrategyN3M5(allc_b));
}

StrategyN3M5 StrategyN3M5::AllD() {
  std::bitset<N> alld_b;
  for (size_t i = 0; i < N; i++) { alld_b.set(i); }
  return std::move(StrategyN3M5(alld_b));
}

StrategyN3M5 StrategyN3M5::TFT() {
  std::bitset<N> tft_b;
  for (size_t i = 0; i < N; i++) {
    const size_t mask_b0 = 1ull << 5ul, mask_c0 = 1ull << 10ul;
    if(i & mask_b0 || i & mask_c0) { tft_b.set(i); }
  }
  return std::move(StrategyN3M5(tft_b));
}

StrategyN3M5 StrategyN3M5::WSLS() {
  std::bitset<N> wsls_b;
  for (size_t i = 0; i < N; i++) {
    const size_t mask_a0 = 1ull << 0ul, mask_b0 = 1ull << 5ul, mask_c0 = 1ull << 10ul;
    const size_t mask = mask_a0 | mask_b0 | mask_c0;
    if ((i & mask) == mask || (i & mask) == 0ull) {} // last action profile is ccc or ddd => C
    else { wsls_b.set(i); }
  }
  return std::move(StrategyN3M5(wsls_b));
}

StrategyN3M5 StrategyN3M5::AON5() {
  std::bitset<N> aon5_b;
  for (size_t i = 0; i < N; i++) {
    const size_t mask_a0 = 0b11111ull << 0ul, mask_b0 = 0b11111ull << 5ul, mask_c0 = 0b11111ull << 10ul;
    const size_t ah = i & mask_a0, bh = (i & mask_b0) >> 5ul, ch = (i & mask_c0) >> 10ul;
    if (ah == bh && ah == ch) { aon5_b.reset(i); }
    else { aon5_b.set(i); }
  }
  return std::move(StrategyN3M5(aon5_b));
}

StrategyN3M5 StrategyN3M5::CAPRI3() {
  typedef std::bitset<15> B;
  auto capri_action_at = [](size_t i)->Action {
    const B I(i);
    std::string s = I.to_string('c', 'D');
    std::reverse(s.begin(), s.end());
    const std::string Istr = s.substr(0,5) + '-' + s.substr(5,5) + '-' + s.substr(10,5);

    const B oldest = 0b10000'10000'10000ul, latest = 0b00001'00001'00001ul;
    const B latest2 = (latest << 1ul) | latest;
    const B a_mask = 0b00000'00000'11111ul;
    const B b_mask = a_mask << 5ul, c_mask = a_mask << 10ul;
    const size_t na = (I & a_mask).count(), nb = (I & b_mask).count(), nc = (I & c_mask).count();

    size_t last_ccc = 5;
    for (size_t t = 0; t < 5; t++) {
      if ((I & (latest<<t)) == B(0ul)) {  // CCC is found at t step before
        last_ccc = t;
        break;
      }
    }

    // C: cooperate if the mutual cooperation is formed at last two rounds
    if ((I & latest2) == 0ul) {
      return C;
    }
      // C0: cooperate if the last action profile is CCC & relative payoff profile is equal
    else if ((I & latest) == 0ul) {
      if (na == nb && na == nc) {
        return C;
      }
    }
    else if (last_ccc > 0 && last_ccc < 5) {
      B mask = latest;
      for (size_t t = 0; t < last_ccc; t++) { mask = ((mask << 1ul) | latest); }
      // A: Accept punishment by prescribing *C* if all your relative payoffs are at least zero.
      size_t pa = (I & mask & a_mask).count();
      size_t pb = (I & mask & b_mask).count();
      size_t pc = (I & mask & c_mask).count();
      if (pa >= pb && pa >= pc) {
        return C;
      }
        // P: Punish by *D* if any of your relative payoffs is negative.
      else {
        return D;
      }
    }

    // R: grab the chance to recover
    if (I == 0b11111'11110'11110 || I == 0b11110'11111'11110 || I == 0b11110'11110'11111) {
      // R: If payoff profile is (+1,+1,-1), prescribe *C*.
      return C;
    }
    if (I == 0b11110'11100'11100 || I == 0b11100'11110'11100 || I == 0b11100'11100'11110) {
      return C;
    }
    // In all other cases, *D*
    return D;
  };

  std::bitset<N> capri3_b;
  for (size_t i = 0; i < N; i++) {
    if (capri_action_at(i) == D) { capri3_b.set(i); }
  }
  return std::move(StrategyN3M5(capri3_b));
}

StrategyN3M5 StrategyN3M5::FUSS_m3() {
  // in m3 results, state i is denoted by a2a1a0_b2b1b0_c2c1c0 in binary
  const std::string m3 = "cdcdcdcdddcdddddcccdcdcdddddddddcdcdcdcdddddddddcdcdcdcdddddddddc"
                         "cddccddcccdcccddcdddcddddddddddccddccddcccdcccddcdddcdddddddddddc"
                         "cddccdcccdccddcccccdccddccddccdccddccdccddccddcdcccdccddccddccccd"
                         "dccddcccdcdcddcddccddddddddcdcccdccddcdcdcdcddcdcdcddddddddddcdcd"
                         "cdcdddddddddcdcdcdcdddddddddcdcdcdcdddddddddcdcdcdcdddddddddccddc"
                         "cddcccdcccddccddcddddddddddccddccddcccdcccddcdddcdddddddddddccddc"
                         "cdccddccddcdcccdccddccddccdccddccdccddccddcdcccdccddccddccccddccd"
                         "dcdcdcdcddcdddcddddddddddccddccddcdcdcdcddcdddcdddddddddd";
  std::bitset<N> fuss_b;
  for (size_t i = 0; i < N; i++) {
    const size_t a_hist = (i & 0b111ul), b_hist = (i&(0b111ul << 5ul)) >> 5ul, c_hist = (i & (7ul << 10ul)) >> 10ul;
    const size_t m3_idx = (a_hist << 6ul) | (b_hist << 3ul) | (c_hist);
    if (m3.at(m3_idx) == 'c') { fuss_b.reset(i); }
    else if (m3.at(m3_idx) == 'd') { fuss_b.set(i); }
    else { throw std::runtime_error("invalid input format"); }
  }
  return std::move(StrategyN3M5(fuss_b));
}
