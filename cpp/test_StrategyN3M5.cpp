#include <iostream>
#include <cassert>
#include "StrategyN3M5.hpp"

#define myassert(x) do {                              \
if (!(x)) {                                           \
  printf("Assertion failed: %s, file %s, line %d\n"   \
         , #x, __FILE__, __LINE__);                   \
  exit(1);                                            \
  }                                                   \
} while (0)

void test_State() {
  StateN3M5 s("dcdcd_ddccd_ccccc");
  myassert(s.ha[0] == D && s.ha[1] == C && s.ha[2] == D && s.ha[3] == C && s.ha[4] == D);
  myassert(s.hb[0] == D && s.hb[1] == D && s.hb[2] == C && s.hb[3] == C && s.hb[4] == D);
  myassert(s.hc[0] == C && s.hc[1] == C && s.hc[2] == C && s.hc[3] == C && s.hc[4] == C);

  myassert(s.ToString() == "dcdcd_ddccd_ccccc");

  myassert(s.NextState(D, C, D) == StateN3M5("ddcdc_cddcc_dcccc"));

  myassert(s.RelativePayoff(true) == 0);
  myassert(s.RelativePayoff(false) == 1);

  myassert(s.ID() == 0b000001001110101);

  myassert(StateN3M5(0b111110101011100) == StateN3M5("ccddd_cdcdc_ddddd"));

  myassert(StateN3M5("ddddd_dcdcd_ccccc").StateFromB() == StateN3M5("dcdcd_ccccc_ddddd"));
  myassert(StateN3M5("ddddd_dcdcd_ccccc").StateFromC() == StateN3M5("ccccc_ddddd_dcdcd"));

  auto noised = StateN3M5("ddddd_dcdcd_ccccc").NoisedStates();
  myassert(noised.size() == 3);
  myassert(std::find(noised.cbegin(), noised.cend(), StateN3M5("cdddd_dcdcd_ccccc")) != noised.end());
  myassert(std::find(noised.cbegin(), noised.cend(), StateN3M5("ddddd_ccdcd_ccccc")) != noised.end());
  myassert(std::find(noised.cbegin(), noised.cend(), StateN3M5("ddddd_dcdcd_dcccc")) != noised.end());

  auto prev = StateN3M5("ddddd_dcdcd_ccccc").PossiblePrevStates();
  myassert(prev.size() == 8);
  std::set<std::string> prev_set;
  for (const auto &st: prev) { prev_set.insert(st.ToString()); }
  std::set<std::string> expected = {"ddddc_cdcdc_ccccc", "ddddc_cdcdc_ccccd", "ddddc_cdcdd_ccccc", "ddddc_cdcdd_ccccd",
                                    "ddddd_cdcdc_ccccc", "ddddd_cdcdc_ccccd", "ddddd_cdcdd_ccccc", "ddddd_cdcdd_ccccd"};
  myassert(prev_set == expected);

  auto s1 = StateN3M5("ddddd_dcdcd_ccccc");
  myassert(s1.NumDiffInT1(StateN3M5("ddddd_dcdcd_ccccc")) == 0);
  myassert(s1.NumDiffInT1(StateN3M5("cdddd_dcdcd_ccccc")) == 1);
  myassert(s1.NumDiffInT1(StateN3M5("cdddd_ccdcd_ccccc")) == 2);
  myassert(s1.NumDiffInT1(StateN3M5("cdddd_ccdcd_dcccc")) == 3);
  myassert(s1.NumDiffInT1(StateN3M5("ddddd_dcdcd_ccccd")) == -1);
}


void test_AllC() {
  std::bitset<32768> allc_b(0ull);
  StrategyN3M5 allc(allc_b);

  std::string s = allc.ToString();
  myassert(s.size() == 32768);
  bool b = true;
  for (char c: s) { if(c != '0') { b = false; break; } }
  myassert(b);

  myassert(allc.ActionAt(StateN3M5(0)) == C);
  myassert(allc.ActionAt(StateN3M5(999)) == C);
  myassert(allc.ActionAt(StateN3M5(9999)) == C);

  myassert(allc.IsEfficientTopo());
  myassert(allc.IsEfficient());

  UnionFind uf = allc.MinimizeDFA(false);
  myassert(uf.to_map().size() == 1);
  myassert(allc.IsDefensibleDFA());
}


void test_AllD() {
  const size_t N = StrategyN3M5::N;
  std::bitset<N> alld_b;
  for (size_t i = 0; i < N; i++) { alld_b.set(i); }
  StrategyN3M5 alld(alld_b);

  myassert(alld.IsDefensibleDFA() == true);
  myassert(alld.IsEfficientTopo() == false);
  myassert(alld.IsEfficient() == false);
  auto dests = alld.DestsOfITG();
  for (int i: dests) { myassert(i == 32767); } // all goes to dddddd

  auto stat = alld.StationaryState(0.00001);
  for (int i = 0; i < 32767; i++) { myassert(stat[i] < 0.01); }
  myassert(stat[32767] > 0.99);

  myassert(alld.IsDistinguishable() == true);
  myassert(alld.IsDistinguishableTopo() == true);

  const auto simp_automaton = alld.MinimizeDFA(false).to_map();
  myassert(simp_automaton.size() == 1);
  const auto full_automaton = alld.MinimizeDFA(true).to_map();
  myassert(full_automaton.size() == 1);
}

void test_TFT() {
  const size_t N = StrategyN3M5::N;
  std::bitset<N> tft_b;
  for (size_t i = 0; i < N; i++) {
    const size_t mask_b0 = 1ull << 5ul, mask_c0 = 1ull << 10ul;
    if(i & mask_b0 || i & mask_c0) { tft_b.set(i); }
  }

  StrategyN3M5 tft(tft_b);
  myassert(tft.IsDefensibleDFA() == true);
  myassert(tft.IsEfficient() == false);
  myassert(tft.IsEfficientTopo() == false);
  auto dests = tft.DestsOfITG();
  for (int i: dests) {
    myassert(i == 0 || i == 32767 );
  } // all goes to either cccccc, dddddd, cdcdcd

  auto assert_near_equal = [](double x, double expected) {
    if (std::abs(x-expected) > 0.01 ) {
      std::cerr << "x: " << x << ", expected: " << expected << std::endl;
      myassert(false);
    }
  };
  auto stat = tft.StationaryState(0.00001);
  assert_near_equal(stat[0], 0.0);
  assert_near_equal(stat[32767], 1.0);

  myassert(tft.IsDistinguishable() == false);
  myassert(tft.IsDistinguishableTopo() == false);

  const auto simp_automaton = tft.MinimizeDFA(false).to_map();
  myassert(simp_automaton.size() == 2);
  myassert(simp_automaton.at(0).size() == 8192);
  myassert(simp_automaton.at(32).size() == 32768-8192);
  const auto full_automaton = tft.MinimizeDFA(true).to_map();
  myassert(full_automaton.size() == 2);
  myassert(full_automaton.at(0).size() == 8192);
  myassert(full_automaton.at(32).size() == 32768-8192);
}

void test_WSLS() {
  const size_t N = StrategyN3M5::N;
  std::bitset<N> wsls_b;
  for (size_t i = 0; i < N; i++) {
    const size_t mask_a0 = 1ull << 0ul, mask_b0 = 1ull << 5ul, mask_c0 = 1ull << 10ul;
    const size_t mask = mask_a0 | mask_b0 | mask_c0;
    if ((i & mask) == mask || (i & mask) == 0ull) {} // last action profile is ccc or ddd => C
    else { wsls_b.set(i); }
  }
  StrategyN3M5 wsls(wsls_b);
  myassert(wsls.IsDefensibleDFA() == false);
  myassert(wsls.IsEfficient() == true);
  myassert(wsls.IsEfficientTopo() == true);
  auto dests = wsls.DestsOfITG();
  for (int i: dests) { myassert(i == 0); } // all goes to cccccc

  auto stat = wsls.StationaryState(0.0001);
  myassert(stat[0] > 0.99);

  myassert(wsls.IsDistinguishable() == true);
  myassert(wsls.IsDistinguishableTopo() == true);

  const auto simp_automaton = wsls.MinimizeDFA(false).to_map();
  myassert(simp_automaton.size() == 2);
  myassert(simp_automaton.at(0).size() == 8192);
  myassert(simp_automaton.at(1).size() == 32768-8192);
  const auto full_automaton = wsls.MinimizeDFA(true).to_map();
  myassert(full_automaton.size() == 2);
  myassert(full_automaton.at(0).size() == 8192);
  myassert(full_automaton.at(1).size() == 32768-8192);
}


void test_m3_FUSS() {
  // in m3 results, state i is denoted by a2a1a0_b2b1b0_c2c1c0 in binary
  const std::string m3 = "cdcdcdcdddcdddddcccdcdcdddddddddcdcdcdcdddddddddcdcdcdcdddddddddc"
                         "cddccddcccdcccddcdddcddddddddddccddccddcccdcccddcdddcdddddddddddc"
                         "cddccdcccdccddcccccdccddccddccdccddccdccddccddcdcccdccddccddccccd"
                         "dccddcccdcdcddcddccddddddddcdcccdccddcdcdcdcddcdcdcddddddddddcdcd"
                         "cdcdddddddddcdcdcdcdddddddddcdcdcdcdddddddddcdcdcdcdddddddddccddc"
                         "cddcccdcccddccddcddddddddddccddccddcccdcccddcdddcdddddddddddccddc"
                         "cdccddccddcdcccdccddccddccdccddccdccddccddcdcccdccddccddccccddccd"
                         "dcdcdcdcddcdddcddddddddddccddccddcdcdcdcddcdddcdddddddddd";
  const size_t N = StrategyN3M5::N;
  std::bitset<N> fuss_b;
  for (size_t i = 0; i < N; i++) {
    const size_t a_hist = (i & 7ul), b_hist = (i&(7ul << 5ul)) >> 5ul, c_hist = (i & (7ul << 10ul)) >> 10ul;
    const size_t m3_idx = (a_hist << 6ul) | (b_hist << 3ul) | (c_hist);
    if (m3.at(m3_idx) == 'c') { fuss_b.reset(i); }
    else if (m3.at(m3_idx) == 'd') { fuss_b.set(i); }
    else { throw std::runtime_error("invalid input format"); }
  }

  StrategyN3M5 fuss(fuss_b);

  const auto simp_automaton = fuss.MinimizeDFA(false).to_map();
  std::cerr << "autom_size: " << simp_automaton.size() << std::endl;
  for (const auto &kv: simp_automaton) {
    std::cerr << kv.first << " => " << kv.second.size() << " [\n  ";
    for (const auto &x: kv.second) {
      std::cerr << x << ", ";
    }
    std::cerr << "]," << std::endl;
  }
  myassert(fuss.IsDefensibleDFA() == true);
  myassert(fuss.IsEfficient() == true);
  myassert(fuss.IsEfficientTopo() == true);

  auto stat = fuss.StationaryState(0.00001);
  myassert(stat[0] > 0.99);
  std::cerr << "stationary state" << std::endl;
  for (size_t i = 0; i < stat.size(); i++) {
    if (stat[i] > 0.05) { std::cerr << StateN3M5(i).ToString() << " : " << stat[i] << std::endl; }
  }

  myassert(fuss.IsDistinguishable() == true);
  myassert(fuss.IsDistinguishableTopo() == true);

  myassert(simp_automaton.size() == 12);
}

/*
void test_TFTATFT() {
  // 0  *cc*cc : c , 16 *dc*cc : c
  // 1  *cc*cd : d , 17 *dc*cd : d
  // 2  *cc*dc : c , 18 *dc*dc : c
  // 3  *cc*dd : d , 19 *dc*dd : c
  // 8  *cd*cc : d , 24 *dd*cc : d
  // 9  *cd*cd : c , 25 *dd*cd : c
  // 10 *cd*dc : c , 26 *dd*dc : c
  // 11 *cd*dd : d , 27 *dd*dd : d
  std::map<int, Action> m = {
      {0, C}, {1, D}, {2, C}, {3, D},
      {8, D}, {9, C}, {10, C}, {11, D},
      {16, C}, {17, D}, {18, C}, {19, C},
      {24, D}, {25, C}, {26, C}, {27, D}
  };
  StrategyN2M3 tft_atft("cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc");
  int mask = 27; // 011011
  for (int i = 0; i < 64; i++) {
    int masked = i & mask;
    tft_atft.SetAction(i, m.at(masked));
  }

  myassert(tft_atft.ToString() == "cdcdcdcddccddccdcdcccdccdccddccdcdcdcdcddccddccdcdcccdccdccddccd");

  myassert(tft_atft.IsDefensible());
  myassert(tft_atft.IsDefensibleDFA());
  myassert(tft_atft.IsEfficient());

  myassert(tft_atft.IsDistinguishable());
  myassert(tft_atft.IsDistinguishableTopo());

  const auto simp_automaton = tft_atft.MinimizeDFA(false).to_map();
  myassert(simp_automaton.size() == 4);
  myassert(simp_automaton.at(0).size() == 28);  // TFT-c
  myassert(simp_automaton.at(1).size() == 20);  // TFT-d
  myassert(simp_automaton.at(8) == std::set<size_t>({8, 12, 40, 44, 24, 28, 56, 60}));  // ATFT-d
  myassert(simp_automaton.at(9) == std::set<size_t>({9, 13, 41, 45, 25, 29, 57, 61}));  // ATFT-c

  const auto full_auto = tft_atft.MinimizeDFA(true).to_map();
  myassert(full_auto.size() == 6);
  myassert(full_auto.at(0).size() == 24);  // TFT-c
  myassert(full_auto.at(1).size() == 12);  // TFT-d
  myassert(full_auto.at(8) == std::set<size_t>({8, 12, 40, 44, 24, 28, 56, 60}));  // ATFT-d
  myassert(full_auto.at(9) == std::set<size_t>({9, 13, 41, 45, 25, 29, 57, 61}));  // ATFT-c
  myassert(full_auto.at(19) == std::set<size_t>({19, 23, 51, 55}));  // TFT-c-2
  myassert(full_auto.at(11) == std::set<size_t>({11, 15, 43, 47, 27, 31, 59, 63}));  // TFT-d-2
}

void test_CAPRI() {
  // action table of CAPRI
  // A's history : B's history (ccc,ccd,cdc,cdd,dcc,dcd,ddc,ddd)
  // ccc : cddd cddd
  // ccd : cdcd dddd
  // cdc : dcdd cddd
  // cdd : dddd dddd
  // dcc : cdcd cdcd
  // dcd : dddd dddd
  // ddc : dddd cdcc
  // ddd : dddd ddcd
  StrategyN2M3 capri("cdddcdddcdcddddddcddcdddddddddddcdcdcdcdddddddddddddcdccddddddcd");

  myassert(capri.IsDefensible());
  myassert(capri.IsDefensibleDFA());
  myassert(capri.IsEfficient());
  myassert(capri.IsEfficientTopo());

  myassert(capri.IsDistinguishable());
  myassert(capri.IsDistinguishableTopo());

  { // distinguishable against WSLS
    StrategyN2M3 wsls("cdcdcdcddcdcdcdccdcdcdcddcdcdcdccdcdcdcddcdcdcdccdcdcdcddcdcdcdc");
    const auto stat = capri.StationaryState(1.0e-6, &wsls);
    myassert(std::abs(stat[61] - 0.5) < 0.01);  // ddd,dcd ~ 0.5
    myassert(std::abs(stat[58] - 0.5) < 0.01);  // ddd,cdc ~ 0.5
    myassert(stat[0] < 0.01);  // ccc,ccc ~ 0.5
  }

  const auto simp_auto = capri.MinimizeDFA(false).to_map();
  myassert(simp_auto.size() == 7);
  const auto full_auto = capri.MinimizeDFA(true).to_map();
  myassert(full_auto.size() == 14);
}

void test_EfficiencyDefensible() {
  StrategyN2M3 s1("cdddddddcdcdddcddcddcdddddccdcddddcdcdddddddddcdddddddcddddcdddd");  // efficient and defensible
  myassert(s1.IsEfficient());
  myassert(s1.IsDefensible());
  myassert(s1.IsDefensibleDFA());
  // auto stat = s1.StationaryState(0.0001);
  // for(int i=0; i<64; i++) { myassert(stat[i] < 0.01); }
}
 */

int main() {
  std::cout << "Testing StrategyN3M5 class" << std::endl;

  // test_State();
  // test_AllC();
  // test_AllD();
  // test_TFT();
  // test_WSLS();
  test_m3_FUSS();
  // test_EfficiencyDefensible();
  // test_TFTATFT();
  // test_CAPRI();
  return 0;
}

