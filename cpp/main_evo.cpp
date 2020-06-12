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
#include "StrategyN3M5.hpp"

// memory-1 species in strategy space discretized with `D`
// memory-1 strategy is characterized by tuple (p_{c0}, p_{c1}, p_{c2}, p_{d0}, p_{d1}, p_{d2}),
// where each of which denotes cooperation probability conditioned by the last Alice's action and the other players' defectors.
// Each of which can take discrete values [0,1/D,2/D,...,D/D].
// ID of a species is given by an integer with base-(D+1).
// ID can take values [0, 11^6-1 = 1771560]
// ID=0 : AllD, ID=ID_MAX-1 : AllC
const size_t DIS = 4ul;
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
};

const size_t Mem1Species::N;

class Ecosystem {
 public:
  Ecosystem(uint64_t seed) : resident(0), rnd(seed) {};
  Mem1Species resident;
  std::mt19937_64 rnd;
  void UpdateResident(double benefit, double cost, uint64_t N, double sigma, double e) {
    std::uniform_real_distribution<double> uni(0.0, ID_MAX);
    Mem1Species mutant = Mem1Species( uni(rnd) );
    double rho = FixationProb(benefit, cost, N, sigma, e, mutant);
    std::uniform_real_distribution<double> uni1(0.0, 1.0);
    if( uni1(rnd) < rho ) {
      resident = mutant;
    }
  }

  double FixationProb(double benefit, double cost, uint64_t N, double sigma, double e, Mem1Species& mutant) {
    // rho_inv = \sum_{i=0}^{N-1} exp(sigma[S]),
    // where S is defined as
    // S =  i/6(i^2−3iN+6i+3N^2−12N+11)s_{yyy}
    //     −i/6(i+1)(2i−3N+4)s_{yyx}
    //     +i/6(i−1)(i+1)s_{yxx}
    //     −i/6(i^2−3i(N−1)+3N^2−6N+2)s_{xyy}
    //     +i/6(i−1)(2i−3N+2))s_{xxy}
    //     −i/6(i^2−3i+2)s_{xxx}

    auto calc_payoff_from_ss = [benefit,cost](const std::array<double,8>& ss) {
      const double s =
            ss[0] * (3*benefit-cost)    // ccc
          + ss[1] * (2*benefit)         // dcc
          + ss[2] * (2*benefit-cost)    // cdc
          + ss[3] * (1*benefit)         // ddc
          + ss[4] * (2*benefit-cost)    // ccd
          + ss[5] * (1* benefit)        // dcd
          + ss[6] * (1*benefit-cost)    // cdd
          + ss[7] * 0.0;                // ddd
      return s;
    };

    const auto xxx = mutant.StationaryState(mutant, mutant, e);
    const double s_xxx = calc_payoff_from_ss(xxx);
    const auto xxy = mutant.StationaryState(mutant, resident, e);
    const double s_xxy = calc_payoff_from_ss(xxy);
    const auto xyy = mutant.StationaryState(resident, resident, e);
    const double s_xyy = calc_payoff_from_ss(xyy);
    const auto yyy = resident.StationaryState(resident, resident, e);
    const double s_yyy = calc_payoff_from_ss(yyy);
    const auto yyx = resident.StationaryState(resident, mutant, e);
    const double s_yyx = calc_payoff_from_ss(yyx);
    const auto yxx = resident.StationaryState(mutant, mutant, e);
    const double s_yxx = calc_payoff_from_ss(yxx);

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

void PrintAbundance(uint64_t tmax, const std::vector<size_t> &histo, const std::string &fname) {
  typedef std::__1::pair<size_t, size_t> SS;
  auto compare_by_val = [](const SS &a, const SS &b) { return (a.second < b.second); };
  std::__1::multiset<std::__1::pair<size_t, size_t>, decltype(compare_by_val)> sorted_histo(compare_by_val);
  for (size_t i = 0; i < histo.size(); i++) {
    if (histo[i] > 0) {
      sorted_histo.insert( std::make_pair(i, histo[i]) );
    }
  }
  std::ofstream fout(fname);
  for (auto it = sorted_histo.rbegin(); it != sorted_histo.rend(); it++) {
    fout << Mem1Species(it->first).ToString() << ' ' << (double)it->second / tmax * ID_MAX << std::endl;
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

  Ecosystem eco(seed);
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

  std::vector<size_t> histo(ID_MAX, 0ul);
  double coop_rate_sum = 0.0;
  for(uint64_t t = 0; t < tmax; t++) {
    eco.UpdateResident(benefit, cost, N, sigma, e);
    size_t res = eco.resident.id;
    histo[res] += 1;

    auto s = eco.resident.StationaryState(eco.resident, eco.resident, e);
    double coop_rate = s[0] * 1.0 + (s[1]+s[2]+s[4]) * (2.0/3.0) + (s[3]+s[5]+s[6]) * (1.0/3.0);
    coop_rate_sum += coop_rate;

    if ( t % t_int == t_int - 1) {
      std::cerr << t << ' ' << coop_rate << std::endl;
    }
    if ( t % t_measure == t_measure -1 ) {
      tout << t << ' ' << coop_rate << std::endl;
    }
  }

  PrintAbundance(tmax, histo, "abundance.dat");

  std::ofstream jout("_output.json");
  double c0 = 0.0, c1 = 0.0, c2 = 0.0, d0 = 0.0, d1 = 0.0, d2 = 0.0;
  for (size_t i = 0; i < histo.size(); i++) {
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

  return 0;
}
