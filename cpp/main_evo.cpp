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
  Mem1Species(size_t _id, size_t _DIS) : id(_id), DIS(_DIS), prob(_id, _DIS) {
    assert(id <= N_M1_Species(DIS));
  }
  size_t id;
  size_t DIS;
  Cprobs prob;

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

  std::vector<double> StationaryState(const Mem1Species& Bstr, const Mem1Species& Cstr, double error = 0.0) const {
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
    std::vector<double> ans(8, 0.0);
    for(int i=0; i<ans.size(); i++) {
      ans[i] = x(i);
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
  Species(size_t ID, size_t DIS) : m1(Mem1Species(0, DIS)), m5(std::bitset<StrategyN3M5::N>()){
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
      m5 = StrategyN3M5::AON(3);
      name = "AON3";
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
  bool is_m1;
  Mem1Species m1;
  StrategyN3M5 m5;
  std::string name;
  std::string ToString() const { return name; }
  std::vector<double> StationaryState(const Species &sb, const Species &sc, double error) const {
    if (is_m1 && sb.is_m1 && sc.is_m1) {
      return m1.StationaryState(sb.m1, sc.m1, error);
    }
    std::cerr << "calculating stationary state: " << name << ' ' << sb.name << ' ' << sc.name << std::endl;

    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletVec;

    for (size_t j = 0; j < StrategyN3M5::N; j++) {
      // calculate transition probability from j to i
      const StateN3M5 sj(j);
      double c_a = _CooperationProb(sj);
      double c_b = sb._CooperationProb(sj.StateFromB());
      double c_c = sc._CooperationProb(sj.StateFromC());
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
    // std::cerr << "  transition matrix has been created" << std::endl;

    Eigen::VectorXd b = Eigen::VectorXd::Zero(S);
    b(S-1) = 1.0;

    Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double> > solver;
    solver.compute(A);
    Eigen::VectorXd x = solver.solve(b);

    // std::cerr << "#iterations:     " << solver.iterations() << std::endl;
    // std::cerr << "estimated error: " << solver.error() << std::endl;

    std::vector<double> ans(S, 0.0);
    for (int i = 0; i < S; i++) { ans[i] = x[i]; }
    return ans;
  }
 private:
  double _CooperationProb(const StateN3M5 &s) const {
    if (is_m1) { return m1.CooperationProb(s); }
    else { return m5.ActionAt(s) == C ? 1.0 : 0.0; }
  }

 public:
  static std::vector<Species> Memory1Species(size_t discrete_level) {
    std::vector<Species> ans;
    size_t n = Mem1Species::N_M1_Species(discrete_level);
    for (size_t i = 0; i < n; i++) {
      ans.emplace_back(i, discrete_level);
    }
    return std::move(ans);
  }
  static std::vector<Species> ReactiveMem1Species(size_t discrete_level) {
    std::vector<Species> ans;
    size_t B = discrete_level + 1;
    size_t n = Mem1Species::N_M1_Species(discrete_level);
    size_t B3 = B * B * B;
    for (size_t i = 0; i < n/B3; i++) {
      size_t id = i * B3 + i;
      ans.emplace_back(id, discrete_level);
    }
    return std::move(ans);
  }
};


class Ecosystem {
 public:
  // Ecosystem(size_t discrete_level, bool with_capri) : DIS(discrete_level) {
  //   N_M1 = Mem1Species::N_M1_Species(DIS);
  //   if (with_capri) { N_SPECIES = N_M1 + 3; } // we add these three species CAPRI3, AON5, FUSS
  //   else { N_SPECIES = N_M1; }
  // };
  Ecosystem(const std::vector<Species> &species_pool, double error) :pool(species_pool), N_SPECIES(species_pool.size()), e(error) {
    CalculateSSCache();
  };
  size_t N_SPECIES;
  std::vector<Species> pool;
  const double e;
  typedef std::vector<double> ss_cache_t;
  std::vector<std::vector<ss_cache_t> > ss_cache;
  // ss_cache[i][i] stores the stationary state when PG game is played by (i,i,i)
  // ss_cache[i][j] stores the stationary state when PG game is played by (i,i,j)

  void CalculateSSCache() {
    ss_cache.resize(N_SPECIES);
    for (size_t i = 0; i < N_SPECIES; i++) {
      ss_cache[i].resize(N_SPECIES);
    }

#pragma omp parallel for schedule(dynamic,1)
    for (size_t I=0; I < N_SPECIES * N_SPECIES; I++) {
      size_t i = I / N_SPECIES;
      size_t j = I % N_SPECIES;
      ss_cache[i][j] = pool[i].StationaryState(pool[i], pool[j], e);
    }
  }

  // payoff of species i and j when the game is played by (i,i,j)
  std::array<double,2> PayoffVersus(size_t i, size_t j, double benefit, double cost) const {
    std::array<double, 3> ans = Payoffs(ss_cache[i][j], benefit, cost);
    assert(ans[0] == ans[1]);
    return {ans[0], ans[2]};
  }

  std::array<double,3> Payoffs(const ss_cache_t &ss, double benefit, double cost) const {
    std::array<double, 3> ans = {0.0, 0.0, 0.0};
    if (ss.size() == 8) {
      for (size_t i = 0; i < 8; i++) {
        double pa = 0.0, pb = 0.0, pc = 0.0;
        if ((i & 1ul) == 0) {
          pa -= cost;
          pb += benefit / 2.0;
          pc += benefit / 2.0;
        }
        if ((i & 2ul) == 0) {
          pb -= cost;
          pc += benefit / 2.0;
          pa += benefit / 2.0;
        }
        if ((i & 4ul) == 0) {
          pc -= cost;
          pa += benefit / 2.0;
          pb += benefit / 2.0;
        }
        ans[0] += ss[i] * pa;
        ans[1] += ss[i] * pb;
        ans[2] += ss[i] * pc;
      }
    }
    else {
      assert(ss.size() == StrategyN3M5::N);
      for (size_t i = 0; i < StrategyN3M5::N; i++) {
        StateN3M5 s(i);
        double pa = 0.0, pb = 0.0, pc = 0.0;
        if (!s.ha[0]) {
          pa -= cost;
          pb += benefit / 2.0;
          pc += benefit / 2.0;
        }
        if (!s.hb[0]) {
          pb -= cost;
          pc += benefit / 2.0;
          pa += benefit / 2.0;
        }
        if (!s.hc[0]) {
          pc -= cost;
          pa += benefit / 2.0;
          pb += benefit / 2.0;
        }
        ans[0] += ss[i] * pa;
        ans[1] += ss[i] * pb;
        ans[2] += ss[i] * pc;
      }
    }
    return ans;
  }

  // calculate the equilibrium distribution exactly by linear algebra
  std::vector<double> CalculateEquilibrium(double benefit, double cost, uint64_t N, double sigma) const {
    Eigen::MatrixXd A(N_SPECIES, N_SPECIES);
    #pragma omp parallel for
    for (size_t ii = 0; ii < N_SPECIES * N_SPECIES; ii++) {
      size_t i = ii / N_SPECIES;
      size_t j = ii % N_SPECIES;
      if (i == j) { A(i, j) = 0.0; continue; }
      double p = FixationProb(benefit, cost, N, sigma, i, j);
      // std::cerr << "Fixation prob of mutant (mutant,resident): " << p << " (" << pool[i].ToString() << ", " << pool[j].ToString() << ")" << std::endl;
      A(i, j) = p * (1.0 / N_SPECIES);
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

  double FixationProb(double benefit, double cost, uint64_t N, double sigma, size_t mutant_idx, size_t resident_idx) const {
    // rho_inv = \sum_{i=0}^{N-1} exp(sigma[S]),
    // where S is defined as
    // S =  i/6(i^2−3iN+6i+3N^2−12N+11)s_{yyy}
    //     −i/6(i+1)(2i−3N+4)s_{yyx}
    //     +i/6(i−1)(i+1)s_{yxx}
    //     −i/6(i^2−3i(N−1)+3N^2−6N+2)s_{xyy}
    //     +i/6(i−1)(2i−3N+2))s_{xxy}
    //     −i/6(i^2−3i+2)s_{xxx}
    // std::cerr << "caclulating fixation prob for mutant " << pool[mutant_idx].ToString() << " against resident " << pool[resident_idx].ToString() << std::endl;

    double s_xxx = PayoffVersus(mutant_idx, mutant_idx, benefit, cost)[0];
    double s_yyy = PayoffVersus(resident_idx, resident_idx, benefit, cost)[0];
    auto xxy = PayoffVersus(mutant_idx, resident_idx, benefit, cost);
    double s_xxy = xxy[0];
    double s_yxx = xxy[1];
    auto yyx = PayoffVersus(resident_idx, mutant_idx, benefit, cost);
    double s_yyx = yyx[0];
    double s_xyy = yyx[1];

    // std::cerr << "s_xxx: " << s_xxx << ", s_xxy: " << s_xxy << ", s_xyy: " << s_xyy << std::endl;
    // std::cerr << "s_yyy: " << s_yyy << ", s_yxx: " << s_yxx << ", s_yyx: " << s_yyx << std::endl;

    double num_games = (N-1) * (N-2) / 2.0;
    double rho_inv = 0.0;
    for (int i=0; i < N; i++) {
      double x = sigma / num_games * (i / 6.0) * (
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
      fout << pool[it->first].ToString() << ' ' << it->second << std::endl;
    }
    fout.close();
  }
  std::vector<std::string> SpeciesNames() const {
    std::vector<std::string> ans;
    for(auto s: pool) {
      ans.emplace_back( s.ToString() );
    }
    return std::move(ans);
  }
  double CooperationLevelSpecies(size_t i) const {
    const ss_cache_t &ss = ss_cache[i][i];
    if (ss.size() == 8) {
      double level = 0.0;
      for (size_t s = 0; s < 8; s++) {
        size_t num_c = 3 - std::bitset<3>(s).count();
        level += ss[s] * (num_c / 3.0);
      };
      return level;
    }
    else {
      double level = 0.0;
      for (size_t s = 0; s < StrategyN3M5::N; s++) {
        StateN3M5 state(s);
        size_t num_c = 3ul;
        if (state.ha[0]) num_c -= 1;
        if (state.hb[0]) num_c -= 1;
        if (state.hc[0]) num_c -= 1;
        level += ss[s] * (num_c / 3.0);
      }
      return level;
    }
  }
  double CooperationLevel(const std::vector<double> &eq_rate) const {
    assert(eq_rate.size() == N_SPECIES);
    double ans = 0.0;
    for (size_t i = 0; i < N_SPECIES; i++) {
      double c_lev = CooperationLevelSpecies(i);
      ans += eq_rate[i] * c_lev;
    }
    return ans;
  }
};


int main(int argc, char *argv[]) {
  Eigen::initParallel();
  if( argc != 5 ) {
    std::cerr << "Error : invalid argument" << std::endl;
    std::cerr << "  Usage : " << argv[0] << " <Nmax> <sigma> <error rate> <discrete_level>" << std::endl;
    return 1;
  }

  double cost = 1.0;
  uint64_t Nmax = std::strtoull(argv[1], nullptr,0);
  double sigma = std::strtod(argv[2], nullptr);
  double e = std::strtod(argv[3], nullptr);
  uint64_t discrete_level = std::strtoull(argv[4], nullptr,0);

  /*
  Species aon(65, 1);
  Species capri3(64, 1);
  auto ss = aon.StationaryState(aon, capri3, e);
  for (size_t i = 0; i < ss.size(); i++) {
    if (ss.at(i) > 0.05) {
      std::cerr << StateN3M5(i).ToString() << " : " << ss.at(i) << std::endl;
    }
  }
  return 0;
   */


  std::vector<Species> pool = Species::ReactiveMem1Species(discrete_level);
  // std::vector<Species> pool = Species::Memory1Species(discrete_level);
  // for (size_t i = 0; i < 64; i++) {
  //   pool.emplace_back(i, discrete_level);
  // }
  pool.emplace_back(64, discrete_level);
  // pool.emplace_back(65, discrete_level);
  // pool.emplace_back(66, discrete_level);
  Ecosystem eco(pool, e);

  {
    std::ofstream namout("species.txt");
    for(auto name: eco.SpeciesNames()) {
      namout << name << std::endl;
    }
    namout << std::endl;
  }

  auto SweepOverBeta = [&eco,cost,sigma](size_t N) {
    char fname1[100], fname2[100];
    sprintf(fname1, "abundance_%zu.dat", N);
    std::ofstream eqout(fname1);
    sprintf(fname2, "cooperation_level_%zu.dat", N);
    std::ofstream coout(fname2);
    for (int i = 5; i <= 300; i+=5) {
      double benefit = 1.0 + i / 100.0;
      auto eq = eco.CalculateEquilibrium(benefit, cost, N, sigma);
      eqout << benefit << ' ';
      for (double x: eq) { eqout << x << ' '; }
      eqout << std::endl;
      double c_lev = eco.CooperationLevel(eq);
      coout << benefit << ' ' << c_lev << std::endl;
    }
  };

  for (int N = 3; N <= Nmax; N++) {
    SweepOverBeta(N);
  }

  return 0;
}
