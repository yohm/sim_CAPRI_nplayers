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
const size_t DIS = 2ul;
const size_t ID_MAX = (DIS+1)*(DIS+1)*(DIS+1)*(DIS+1)*(DIS+1)*(DIS+1);
class Cprobs {
 public:
  Cprobs(size_t id) {
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
  static const size_t N = 8; // 2^{nm} = 8
  Mem1Species(size_t _id) : id(_id), prob(_id), _ss_cache_ready(false) {
    assert(id <= ID_MAX);
  }
  size_t id;
  Cprobs prob;
  std::array<double,N> _ss_cache;
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

  std::array<double,N> StationaryState(const Mem1Species& Bstr, const Mem1Species& Cstr, double error = 0.0) {
    if (id == Bstr.id && id == Cstr.id && _ss_cache_ready) { return _ss_cache; }
    Eigen::Matrix<double,N,N> A;

    // state 0: ccc, state 1: ccd (last bit is A's history), ... 7: ddd
    // calculate transition probability from j to i

    for (size_t j = 0; j < N; j++) {
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

    for(int i=0; i<N; i++) {
      A(i,i) -= 1.0;
    }
    for(int i=0; i<N; i++) {
      A(N-1,i) += 1.0;  // normalization condition
    }
    Eigen::VectorXd b(N);
    for(int i=0; i<N; i++) { b(i) = 0.0;}
    b(N-1) = 1.0;
    Eigen::VectorXd x = A.colPivHouseholderQr().solve(b);
    std::array<double,N> ans = {0};
    for(int i=0; i<N; i++) {
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

const size_t Mem1Species::N;

class Species { // either Mem1Species or StrategyN3M5
 public:
  Species(size_t ID) : m1(Mem1Species(0)), m5(std::bitset<StrategyN3M5::N>()) {
    if (ID < ID_MAX) {
      is_m1 = true;
      m1 = Mem1Species(ID);
      m5 = StrategyN3M5(std::bitset<StrategyN3M5::N>());
    }
    else if (ID == ID_MAX) {
      is_m1 = false;
      m1 = Mem1Species(0);
      m5 = StrategyN3M5::CAPRI3();
    }
    else {
      throw std::runtime_error("must not happen");
    }
  };
  Species(const Mem1Species &_m1) : is_m1(true), m1(_m1), m5(StrategyN3M5(std::bitset<StrategyN3M5::N>())) {};
  Species(const StrategyN3M5 &_m5) : is_m1(false), m1(Mem1Species(0)), m5(_m5) {};
  Species & operator=(const Species & rhs) { is_m1 = rhs.is_m1; m1 = rhs.m1; m5 = rhs.m5; return *this; }
  bool is_m1;
  Mem1Species m1;
  StrategyN3M5 m5;
  std::string ToString() const {
    if (is_m1) { return m1.ToString(); }
    else { return std::string("CAPRI3"); }
  }
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

    // Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > solver;
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

};

typedef std::array<double, StrategyN3M5::N> ProbArray;


class Ecosystem {
 public:
  Ecosystem(uint64_t seed, bool with_capri) : resident(Mem1Species(0)), rnd(seed) {
    if (with_capri) { N_SPECIES = ID_MAX + 1; }
    else { N_SPECIES = ID_MAX; }
  };
  Species resident;
  std::mt19937_64 rnd;
  size_t N_SPECIES;
  size_t ResidentID() const {
    if (resident.is_m1) return resident.m1.id;
    else return ID_MAX;
  }
  void UpdateResident(double benefit, double cost, uint64_t N, double sigma, double e) {
    std::uniform_int_distribution<size_t> uni(0, N_SPECIES-1);
    size_t r = uni(rnd);
    Species mutant = (r < ID_MAX) ? Species(Mem1Species(r)) : Species(StrategyN3M5::CAPRI3());

    double rho = FixationProb(benefit, cost, N, sigma, e, mutant);
    std::uniform_real_distribution<double> uni1(0.0, 1.0);
    if( uni1(rnd) < rho ) {
      resident = mutant;
    }
  }

  // calculate the equilibrium distribution exactly by linear algebra
  std::vector<double> CalculateEquilibrium(double benefit, double cost, uint64_t N, double sigma, double e) {
    Eigen::initParallel();
    Eigen::MatrixXd A(N_SPECIES, N_SPECIES);
    #pragma omp parallel for
    for (size_t i = 0; i < ID_MAX; i++) {
      Species si = (i < ID_MAX) ? Species(Mem1Species(i)) : Species(StrategyN3M5::CAPRI3());
      for (size_t j = 0; j < N_SPECIES; j++) {
        if (i == j) continue;
        // calculate the transition probability from j to i
        Species sj = (j < ID_MAX) ? Species(Mem1Species(j)) : Species(StrategyN3M5::CAPRI3());
        if (i == ID_MAX || j == ID_MAX) {
          std::cerr << "calculating rho for (" << i << ", " << j << ")" << std::endl;
        }
        double p = FixationProb(benefit, cost, N, sigma, e, si, sj);
        A(i, j) = p * (1.0 / N_SPECIES);
      }
    }
    for (size_t i = ID_MAX; i < N_SPECIES; i++) {
      Species si = (i < ID_MAX) ? Species(Mem1Species(i)) : Species(StrategyN3M5::CAPRI3());
      #pragma omp parallel for
      for (size_t j = 0; j < N_SPECIES; j++) {
        if (i == j) continue;
        // calculate the transition probability from j to i
        Species sj = (j < ID_MAX) ? Species(Mem1Species(j)) : Species(StrategyN3M5::CAPRI3());
        if (i == ID_MAX || j == ID_MAX) {
          std::cerr << "calculating rho for (" << i << ", " << j << ")" << std::endl;
        }
        double p = FixationProb(benefit, cost, N, sigma, e, si, sj);
        A(i, j) = p * (1.0 / N_SPECIES);
      }
    }
    for (size_t i = 0; i < N_SPECIES; i++) {
      double p_sum = 0.0;
      for (size_t j = 0; j < N_SPECIES; j++) {
        if (i == j) continue;
        p_sum += A(j, i);
      }
      A(i, i) = 1.0 - p_sum; // probability that the state doesn't change
    }

    // subtract Ax = x => (A-I)x = 0
    for (size_t i = 0; i < N_SPECIES; i++) {
      A(i, i) -= 1.0;
    }
    // normalization condition
    for (size_t i = 0; i < N_SPECIES; i++) {
      A(N_SPECIES-1, i) += 1.0;
    }

    Eigen::VectorXd b(N_SPECIES);
    for(int i=0; i<N_SPECIES; i++) { b(i) = 0.0;}
    b(N_SPECIES-1) = 1.0;
    Eigen::VectorXd x = A.householderQr().solve(b);
    std::vector<double> ans(N_SPECIES);
    for(int i=0; i<N_SPECIES; i++) {
      ans[i] = x(i);
    }
    return ans;
  }

  double FixationProb(double benefit, double cost, uint64_t N, double sigma, double e, Species& mutant) {
    return FixationProb(benefit, cost, N, sigma, e, mutant, resident);
  }

  double FixationProb(double benefit, double cost, uint64_t N, double sigma, double e, Species& mutant, Species & res) {
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
    const auto xxy = mutant.Payoffs(mutant, res, benefit, cost, e);
    const double s_xxy = xxy[0];
    const auto xyy = mutant.Payoffs(res, res, benefit, cost, e);
    const double s_xyy = xyy[0];
    const auto yyy = res.Payoffs(res, res, benefit, cost, e);
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
};

void PrintAbundance(const std::vector<double> &histo_d, const std::string &fname) {
  typedef std::pair<size_t, double> SS;
  auto compare_by_val = [](const SS &a, const SS &b) { return (a.second < b.second); };
  std::multiset<SS, decltype(compare_by_val)> sorted_histo(compare_by_val);
  for (size_t i = 0; i < histo_d.size(); i++) {
    if (histo_d[i] > 0.0) {
      sorted_histo.insert( std::make_pair(i, histo_d[i]) );
    }
  }
  std::ofstream fout(fname);
  for (auto it = sorted_histo.rbegin(); it != sorted_histo.rend(); it++) {
    fout << Species(it->first).ToString() << ' ' << it->second << std::endl;
  }
  fout.close();
}

int main(int argc, char *argv[]) {
  if( argc != 9 ) {
    std::cerr << "Error : invalid argument" << std::endl;
    std::cerr << "  Usage : " << argv[0] << " <benefit> <cost> <N> <N_sigma> <error rate> <t_init> <tmax> <rand_seed>" << std::endl;
    return 1;
  }

  double benefit = std::strtod(argv[1], nullptr);
  double cost = std::strtod(argv[2], nullptr);
  uint64_t N = std::strtoull(argv[3], nullptr,0);
  double N_sigma = std::strtod(argv[4], nullptr);
  double sigma = N_sigma / N;
  double e = std::strtod(argv[5], nullptr);
  uint64_t t_init = std::strtoull(argv[6], nullptr,0);
  uint64_t tmax = std::strtoull(argv[7], nullptr,0);
  uint64_t seed = std::strtoull(argv[8], nullptr,0);

  Ecosystem eco(seed, true);
  uint64_t t_int = 10000;
  for (uint64_t t = 0; t < t_init; t++) {
    eco.UpdateResident(benefit, cost, N, sigma, e);
    if( t % t_int == t_int - 1) {
      double c = eco.resident.StationaryState(eco.resident, eco.resident, e)[0];
      std::cerr << t << ' ' << c << std::endl;
    }
  }

  std::ofstream tout("timeseries.dat");
  uint64_t t_measure = tmax / 100000 + 1;

  std::vector<size_t> histo(eco.N_SPECIES, 0ul);
  double coop_rate_sum = 0.0;
  for(uint64_t t = 0; t < tmax; t++) {
    eco.UpdateResident(benefit, cost, N, sigma, e);
    size_t res = eco.ResidentID();
    histo[res] += 1;

    // auto s = eco.resident.StationaryState(eco.resident, eco.resident, e);
    // double coop_rate = s[0] * 1.0 + (s[1]+s[2]+s[4]) * (2.0/3.0) + (s[3]+s[5]+s[6]) * (1.0/3.0);
    // coop_rate_sum += coop_rate;

    if ( t % t_int == t_int - 1) {
      // std::cerr << t << ' ' << coop_rate << std::endl;
    }
    if ( t % t_measure == t_measure -1 ) {
      // tout << t << ' ' << coop_rate << std::endl;
    }
  }

  std::vector<double> histo_d(eco.N_SPECIES);
  for (size_t i = 0; i < histo.size(); i++) {
    histo_d[i] = (double)histo[i] / tmax * eco.N_SPECIES;
  }
  PrintAbundance(histo_d, "abundance.dat");

  std::ofstream jout("_output.json");
  double c0 = 0.0, c1 = 0.0, c2 = 0.0, d0 = 0.0, d1 = 0.0, d2 = 0.0;
  for (size_t i = 0; i < ID_MAX; i++) {
    Mem1Species s(i);
    c0 += histo[i] * s.prob.c0;
    c1 += histo[i] * s.prob.c1;
    c2 += histo[i] * s.prob.c2;
    d0 += histo[i] * s.prob.d0;
    d1 += histo[i] * s.prob.d1;
    d2 += histo[i] * s.prob.d2;
  }
  jout << "{" << std::endl;
  jout << "\"coop_rate\": " << coop_rate_sum / tmax << "," << std::endl;
  jout << "\"c0\": " << c0 / tmax << "," << std::endl;
  jout << "\"c1\": " << c1 / tmax << "," << std::endl;
  jout << "\"c2\": " << c2 / tmax << "," << std::endl;
  jout << "\"d0\": " << d0 / tmax << "," << std::endl;
  jout << "\"d1\": " << d1 / tmax << "," << std::endl;
  jout << "\"d2\": " << d2 / tmax << "\n}" << std::endl;
  jout.close();

  std::cerr << "Calculating equilibrium" << std::endl;
  auto eq = eco.CalculateEquilibrium(benefit, cost, N, sigma, e);
  PrintAbundance(eq, "equilibrium.dat");

  return 0;
}
