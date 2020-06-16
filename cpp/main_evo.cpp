//
// Created by Yohsuke Murase on 2020/06/04.
//

#include <iostream>
#include <vector>
#include <array>
#include <map>
#include <random>
#include <cassert>
#include <fstream>
#include <Eigen/Dense>
#include "StrategyN3M5.hpp"

// memory-1 species in strategy space discretized with `D`
// memory-1 strategy is characterized by tuple (p_{c0}, p_{c1}, p_{c2}, p_{d0}, p_{d1}, p_{d2}),
// where each of which denotes cooperation probability conditioned by the last Alice's action and the other players' defectors.
// Each of which can take discrete values [0,1/D,2/D,...,D/D].
// ID of a species is given by an integer with base-(D+1).
// ID can take values [0, (D+1)^6-1]
// ID=0 : AllD, ID=ID_MAX-1 : AllC
class Cprobs {
 public:
  Cprobs(size_t id, size_t DIS) {
    double dinv = 1.0 / DIS;
    const size_t B = DIS+1;
    c0 = (id % B) * dinv;
    c1 = ((id/B) % B) * dinv;
    c2 = ((id/B/B) % B) * dinv;
    d0 = ((id/B/B/B) % B) * dinv;
    d1 = ((id/B/B/B/B) % B) * dinv;
    d2 = ((id/B/B/B/B/B) % B) * dinv;
  }
  double c0, c1, c2, d0, d1, d2;

};

class Mem1Species {
 public:
  static size_t N_M1_Species(size_t DIS) { return (DIS+1)*(DIS+1)*(DIS+1)*(DIS+1)*(DIS+1)*(DIS+1); }
  Mem1Species(size_t _id, size_t _DIS) : id(_id), DIS(_DIS), prob(_id, _DIS), _ss_cache_ready(false) {
    assert(id <= N_M1_Species(DIS));
  }
  size_t id;
  size_t DIS;
  Cprobs prob;
  std::array<double,8> _ss_cache;
  bool _ss_cache_ready;

  std::string ToString() const {
    const size_t B = DIS+1;
    size_t c0 = (id % B);
    size_t c1 = ((id/B) % B);
    size_t c2 = ((id/B/B) % B);
    size_t d0 = ((id/B/B/B) % B);
    size_t d1 = ((id/B/B/B/B) % B);
    size_t d2 = ((id/B/B/B/B/B) % B);
    std::ostringstream oss;
    oss << c0 << '-' << c1 << '-' << c2 << '-' << d0 << '-' << d1 << '-' << d2;
    return oss.str();
  }

  std::array<double,8> StationaryState(const Mem1Species& Bstr, const Mem1Species& Cstr, double error = 0.0) {
    if (id == Bstr.id && id == Cstr.id && _ss_cache_ready) { return _ss_cache; }
    Eigen::Matrix<double,8,8> A;

    // state 0: ccc, state 1: ccd (last bit is A's history), ... 7: ddd
    // calculate transition probability from j to i

    for (size_t j = 0; j < 8; j++) {
      Action last_a = (j & 1ul) ? D : C;
      Action last_b = ((j>>1ul) & 1ul) ? D : C;
      Action last_c = ((j>>2ul) & 1ul) ? D : C;

      auto cooperation_prob = [](Action a, Action b, Action c, const Mem1Species& str) -> double {
        if (a == C) {
          if (b == C && c == C) { return str.prob.c0; }
          else if (b == D && c == C) { return str.prob.c1; }
          else if (b == C && c == D) { return str.prob.c1; }
          else { return str.prob.c2; }
        }
        else {
          if (b == C && c == C) { return str.prob.d0; }
          else if (b == D && c == C) { return str.prob.d1; }
          else if (b == C && c == D) { return str.prob.d1; }
          else { return str.prob.d2; }
        }
      };
      double c_A = cooperation_prob(last_a, last_b, last_c, *this);
      double c_B = cooperation_prob(last_b, last_c, last_a, Bstr);
      double c_C = cooperation_prob(last_c, last_a, last_b, Cstr);

      c_A = (1.0 - error) * c_A + error * (1.0 - c_A);
      c_B = (1.0 - error) * c_B + error * (1.0 - c_B);
      c_C = (1.0 - error) * c_C + error * (1.0 - c_C);
      A(0,j) = c_A * c_B * c_C;
      A(1,j) = (1.0-c_A) * c_B * c_C;
      A(2,j) = c_A * (1.0-c_B) * c_C;
      A(3,j) = (1.0-c_A) * (1.0-c_B) * c_C;
      A(4,j) = c_A * c_B * (1.0-c_C);
      A(5,j) = (1.0-c_A) * c_B * (1.0-c_C);
      A(6,j) = c_A * (1.0-c_B) * (1.0-c_C);
      A(7,j) = (1.0-c_A) * (1.0-c_B) * (1.0-c_C);
    }

    for(int i=0; i<8; i++) {
      A(i,i) -= 1.0;
    }
    for(int i=0; i<8; i++) {
      A(8-1,i) += 1.0;  // normalization condition
    }
    Eigen::VectorXd b(8);
    for(int i=0; i<8; i++) { b(i) = 0.0;}
    b(8-1) = 1.0;
    Eigen::VectorXd x = A.colPivHouseholderQr().solve(b);
    std::array<double,8> ans = {0};
    for(int i=0; i<ans.size(); i++) {
      ans[i] = x(i);
    }

    if (id == Bstr.id && id == Cstr.id) {
      std::copy(ans.begin(), ans.end(), _ss_cache.begin());
    }

    return ans;
  }

  double CooperationProb(StateN3M5 s) const {
    size_t num_d = 0;
    if (s.hb[0]) num_d++;
    if (s.hc[0]) num_d++;
    if (s.ha[0]) {  // A's last action is d
      if (num_d == 2) { return prob.d2; }
      else if (num_d == 1) { return prob.d1; }
      else { return prob.d0; }  // num_d == 0
    }
    else {
      if (num_d == 2) { return prob.c2; }
      else if (num_d == 1) { return prob.c1; }
      else { return prob.c0; } // num_d == 0
    }
  };
};


class Species { // either Mem1Species or StrategyN3M5
 public:
  Species(size_t ID, size_t DIS) : m1(Mem1Species(0, DIS)), m5(std::bitset<StrategyN3M5::N>()) {
    const size_t N_M1 = Mem1Species::N_M1_Species(DIS);
    if (ID < N_M1) {
      is_m1 = true;
      m1 = Mem1Species(ID, DIS);
      m5 = StrategyN3M5(std::bitset<StrategyN3M5::N>());
      name = m1.ToString();
    }
    else if (ID == N_M1) {
      is_m1 = false;
      m1 = Mem1Species(0, DIS);
      m5 = StrategyN3M5::CAPRI3();
      name = "CAPRI3";
    }
    else if (ID == N_M1 + 1) {
      is_m1 = false;
      m1 = Mem1Species(0, DIS);
      m5 = StrategyN3M5::AON5();
      name = "AON5";
    }
    else if (ID == N_M1 + 2) {
      is_m1 = false;
      m1 = Mem1Species(0, DIS);
      m5 = StrategyN3M5::FUSS_m3();
      name = "FUSS_m3";
    }
    else {
      throw std::runtime_error("must not happen");
    }
  };
  // Species(const Mem1Species &_m1) : is_m1(true), m1(_m1), m5(StrategyN3M5(std::bitset<StrategyN3M5::N>())) {};
  // Species(const StrategyN3M5 &_m5) : is_m1(false), m1(Mem1Species(0, 1)), m5(_m5) {};
  bool is_m1;
  Mem1Species m1;
  StrategyN3M5 m5;
  std::string name;
  std::string ToString() const { return name; }
  double CooperationProb(const StateN3M5 &s) const {
    if (is_m1) { return m1.CooperationProb(s); }
    else { return m5.ActionAt(s) == C ? 1.0 : 0.0; }
  }
  std::array<double, StrategyN3M5::N> StationaryState(const Species &sb, const Species &sc, double error) {
    std::cerr << "calculating stationary state" << std::endl;

    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletVec;

    for (size_t j = 0; j < StrategyN3M5::N; j++) {
      // calculate transition probability from j to i
      const StateN3M5 sj(j);
      double c_a = CooperationProb(sj);
      double c_b = sb.CooperationProb(sj.StateFromB());
      double c_c = sc.CooperationProb(sj.StateFromC());
      // cooperation probability taking noise into account
      c_a = (1.0 - error) * c_a + error * (1.0 - c_a);
      c_b = (1.0 - error) * c_b + error * (1.0 - c_b);
      c_c = (1.0 - error) * c_c + error * (1.0 - c_c);

      for (size_t t = 0; t < 8; t++) {
        Action act_a = (t & 1ul) ? D : C;
        Action act_b = (t & 2ul) ? D : C;
        Action act_c = (t & 4ul) ? D : C;
        size_t i = sj.NextState(act_a, act_b, act_c).ID();
        double p_a = (act_a == C) ? c_a : (1.0-c_a);
        double p_b = (act_b == C) ? c_b : (1.0-c_b);
        double p_c = (act_c == C) ? c_c : (1.0-c_c);
        tripletVec.emplace_back(i, j, p_a * p_b * p_c);
      }
    }

    const size_t S = StrategyN3M5::N;
    Eigen::SparseMatrix<double> A(S, S);
    A.setFromTriplets(tripletVec.cbegin(), tripletVec.cend());

    // subtract unit matrix & normalization condition
    std::vector<T> iVec;
    for (int i = 0; i < S-1; i++) { iVec.emplace_back(i, i, -1.0); }
    for (int i = 0; i < S-1; i++) { iVec.emplace_back(S-1, i, 1.0); }
    Eigen::SparseMatrix<double> I(S, S);
    I.setFromTriplets(iVec.cbegin(), iVec.cend());
    A = A + I;
    std::cerr << "  transition matrix has been created" << std::endl;

    Eigen::VectorXd b = Eigen::VectorXd::Zero(S);
    b(S-1) = 1.0;

    Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double> > solver;
    solver.compute(A);
    Eigen::VectorXd x = solver.solve(b);

    std::cerr << "#iterations:     " << solver.iterations() << std::endl;
    std::cerr << "estimated error: " << solver.error() << std::endl;

    std::array<double, S> ans = {0};
    for (int i = 0; i < S; i++) { ans[i] = x[i]; }
    return ans;
  }

  std::array<double,3> Payoffs(const Species &sb, const Species &sc, double benefit, double cost, double error) {
    std::array<double, 3> ans = {0.0, 0.0, 0.0};
    if (is_m1 && sb.is_m1 && sc.is_m1) {
      auto ss = m1.StationaryState(sb.m1, sc.m1, error);
      for (size_t i = 0; i < 8; i++) {
        size_t num_c = 0;
        double cost_A = 0.0, cost_B = 0.0, cost_C = 0.0;
        if ((i & 1ul) == 0) { num_c += 1; cost_A += cost; }
        if ((i & 2ul) == 0) { num_c += 1; cost_B += cost; }
        if ((i & 4ul) == 0) { num_c += 1; cost_C += cost; }
        ans[0] += ss[i] * (num_c * benefit - cost_A);
        ans[1] += ss[i] * (num_c * benefit - cost_B);
        ans[2] += ss[i] * (num_c * benefit - cost_C);
      };
    }
    else {
      auto ss = StationaryState(sb, sc, error);
      for (size_t i = 0; i < ss.size(); i++) {
        StateN3M5 s(i);
        size_t num_c = 0;
        double cost_A = 0.0, cost_B = 0.0, cost_C = 0.0;
        if (!s.ha[0]) { num_c += 1; cost_A += cost; }
        if (!s.hb[0]) { num_c += 1; cost_B += cost; }
        if (!s.hc[0]) { num_c += 1; cost_C += cost; }
        ans[0] += ss[i] * (num_c * benefit - cost_A);
        ans[1] += ss[i] * (num_c * benefit - cost_B);
        ans[2] += ss[i] * (num_c * benefit - cost_C);
      }
    }
    return ans;
  }

  double CooperationLevel(double error) {
    if (is_m1) {
      std::array<double,8> ss = m1.StationaryState(m1, m1, error);
      double level = 0.0;
      for (size_t i = 0; i < 8; i++) {
        size_t num_c = 3 - std::bitset<3>(i).count();
        level += ss[i] * (num_c / 3.0);
      };
      return level;
    }
    else {
      std::array<double, StrategyN3M5::N> ss = m5.StationaryState();
      double level = 0.0;
      for (size_t i = 0; i < StrategyN3M5::N; i++) {
        StateN3M5 s(i);
        size_t num_c = 3ul;
        if (s.ha[0]) num_c -= 1;
        if (s.hb[0]) num_c -= 1;
        if (s.hc[0]) num_c -= 1;
        level += ss[i] * (num_c / 3.0);
      }
      return level;
    }
  }
};


class Ecosystem {
 public:
  Ecosystem(size_t discrete_level, bool with_capri) : DIS(discrete_level) {
    N_M1 = Mem1Species::N_M1_Species(DIS);
    if (with_capri) { N_SPECIES = N_M1 + 3; } // we add these three species CAPRI3, AON5, FUSS
    else { N_SPECIES = N_M1; }
  };
  size_t DIS;
  size_t N_SPECIES;
  size_t N_M1;
  // calculate the equilibrium distribution exactly by linear algebra
  std::vector<double> CalculateEquilibrium(double benefit, double cost, uint64_t N, double sigma, double e) {
    Eigen::initParallel();
    std::cerr << "DIS, N_SPECIES, N_M1: " << DIS << ", " << N_SPECIES << ", " << N_M1 << std::endl;
    Eigen::MatrixXd A(N_SPECIES, N_SPECIES);
    #pragma omp parallel for
    for (size_t i = 0; i < N_M1; i++) {
      Species si(i, DIS);
      for (size_t j = 0; j < N_SPECIES; j++) {
        if (i == j) { A(i, j) = 0.0; continue; }
        // calculate the transition probability from j to i
        Species sj(j, DIS);
        if (i >= N_M1 || j >= N_M1) { std::cerr << "calculating rho for (" << i << ", " << j << ")" << std::endl; }
        double p = FixationProb(benefit, cost, N, sigma, e, si, sj);
        A(i, j) = p * (1.0 / N_SPECIES);
      }
    }
    for (size_t i = N_M1; i < N_SPECIES; i++) {
      Species si(i, DIS);
      #pragma omp parallel for
      for (size_t j = 0; j < N_SPECIES; j++) {
        if (i == j) { A(i, j) = 0.0; continue; }
        // calculate the transition probability from j to i
        Species sj(j, DIS);
        if (i >= N_M1 || j >= N_M1) { std::cerr << "calculating rho for (" << i << ", " << j << ")" << std::endl; }
        double p = FixationProb(benefit, cost, N, sigma, e, si, sj);
        assert( p >= 0.0 && p <= 1.0 );
        A(i, j) = p * (1.0 / N_SPECIES);
      }
    }
    for (size_t j = 0; j < N_SPECIES; j++) {
      double p_sum = 0.0;
      for (size_t i = 0; i < N_SPECIES; i++) {
        p_sum += A(i, j);
      }
      assert(p_sum <= 1.0);
      A(j, j) = 1.0 - p_sum; // probability that the state doesn't change
    }

    size_t n_row = A.rows();

    // subtract Ax = x => (A-I)x = 0
    for (size_t i = 0; i < A.rows(); i++) {
      A(i, i) -= 1.0;
    }
    // normalization condition
    for (size_t i = 0; i < A.rows(); i++) {
      A(A.rows()-1, i) += 1.0;
    }

    Eigen::VectorXd b(A.rows());
    for(int i=0; i<A.rows()-1; i++) { b(i) = 0.0;}
    b(A.rows()-1) = 1.0;
    Eigen::VectorXd x = A.householderQr().solve(b);
    std::vector<double> ans(A.rows());
    double prob_total = 0.0;
    for(int i=0; i<ans.size(); i++) {
      ans[i] = x(i);
      prob_total += x(i);
      assert(x(i) > -0.000001);
    }
    assert(std::abs(prob_total - 1.0) < 0.00001);
    return ans;
  }

  double FixationProb(double benefit, double cost, uint64_t N, double sigma, double e, Species& mutant, Species & resident) {
    // rho_inv = \sum_{i=0}^{N-1} exp(sigma[S]),
    // where S is defined as
    // S =  i/6(i^2−3iN+6i+3N^2−12N+11)s_{yyy}
    //     −i/6(i+1)(2i−3N+4)s_{yyx}
    //     +i/6(i−1)(i+1)s_{yxx}
    //     −i/6(i^2−3i(N−1)+3N^2−6N+2)s_{xyy}
    //     +i/6(i−1)(2i−3N+2))s_{xxy}
    //     −i/6(i^2−3i+2)s_{xxx}

    const auto xxx = mutant.Payoffs(mutant, mutant, benefit, cost, e);
    const double s_xxx = xxx[0];
    const auto xxy = mutant.Payoffs(mutant, resident, benefit, cost, e);
    const double s_xxy = xxy[0];
    const auto xyy = mutant.Payoffs(resident, resident, benefit, cost, e);
    const double s_xyy = xyy[0];
    const auto yyy = resident.Payoffs(resident, resident, benefit, cost, e);
    const double s_yyy = yyy[0];
    const double s_yyx = xyy[2];
    const double s_yxx = xxy[2];

    double rho_inv = 0.0;
    for (int i=0; i < N; i++) {
      double x = sigma * (i / 6.0) * (
          (i*i - 3.0*i*N + 6.0*i + 3.0*N*N - 12.0*N + 11.0) * s_yyy
          - (i+1.0) * (2.0*i - 3.0*N + 4) * s_yyx
          + (i-1.0) * (i+1.0) * s_yxx
          - (i*i - 3.0*i*(N-1.0) + 3.0*N*N - 6.0*N + 2.0) * s_xyy
          + (i-1.0) * (2.0*i-3.0*N+2.0) * s_xxy
          - (i*i - 3.0*i + 2.0) * s_xxx
      );
      rho_inv += std::exp(x);
    }
    return 1.0 / rho_inv;
  }

  void PrintAbundance(const std::vector<double> &eq_rate, const std::string &fname) const {
    assert(eq_rate.size() == N_SPECIES);
    typedef std::pair<size_t, double> SS;
    auto compare_by_val = [](const SS &a, const SS &b) { return (a.second < b.second); };
    std::multiset<SS, decltype(compare_by_val)> sorted_histo(compare_by_val);
    for (size_t i = 0; i < eq_rate.size(); i++) {
      sorted_histo.insert( std::make_pair(i, eq_rate[i]) );
    }
    std::ofstream fout(fname);
    for (auto it = sorted_histo.rbegin(); it != sorted_histo.rend(); it++) {
      fout << Species(it->first, DIS).ToString() << ' ' << it->second << std::endl;
    }
    fout.close();
  }
  double CooperationLevel(const std::vector<double> &eq_rate, double error) const {
    assert(eq_rate.size() == N_SPECIES);
    double ans = 0.0;
    for (size_t i = 0; i < N_SPECIES; i++) {
      Species s(i, DIS);
      double c_lev = s.CooperationLevel(error);
      ans += eq_rate[i] * c_lev;
    }
    return ans;
  }
};


int main(int argc, char *argv[]) {
  if( argc != 7 ) {
    std::cerr << "Error : invalid argument" << std::endl;
    std::cerr << "  Usage : " << argv[0] << " <benefit> <cost> <N> <N_sigma> <error rate> <discrete_level>" << std::endl;
    return 1;
  }

  double benefit = std::strtod(argv[1], nullptr);
  double cost = std::strtod(argv[2], nullptr);
  uint64_t N = std::strtoull(argv[3], nullptr,0);
  double N_sigma = std::strtod(argv[4], nullptr);
  double sigma = N_sigma / N;
  double e = std::strtod(argv[5], nullptr);
  uint64_t discrete_level = std::strtoull(argv[6], nullptr,0);

  Ecosystem eco(discrete_level, true);


  std::cerr << "Calculating equilibrium" << std::endl;
  auto eq = eco.CalculateEquilibrium(benefit, cost, N, sigma, e);
  eco.PrintAbundance(eq, "equilibrium.dat");
  double c_lev = eco.CooperationLevel(eq, e);
  std::ofstream jout("_output.json");
  jout << "{ \"cooperation_level\": " << c_lev << " }" << std::endl;
  jout.close();

  return 0;
}
