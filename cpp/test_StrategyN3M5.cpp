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
  StrategyN3M5 allc = StrategyN3M5::AllC();

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

  auto simp_auto = allc.MinimizeDFAHopcroft(false);
  myassert(simp_auto.size() == 1);
  auto full_auto = allc.MinimizeDFAHopcroft(true);
  myassert(full_auto.size() == 1);

  myassert(allc.IsDefensibleDFA() == false);
}


void test_AllD() {
  const StrategyN3M5 alld = StrategyN3M5::AllD();

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

  const auto simp_automaton = alld.MinimizeDFAHopcroft(false).to_map();
  myassert(simp_automaton.size() == 1);
  const auto full_automaton = alld.MinimizeDFAHopcroft(true).to_map();
  myassert(full_automaton.size() == 1);
}

void test_TFT() {
  const StrategyN3M5 tft = StrategyN3M5::TFT();
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

  const auto simp_automaton = tft.MinimizeDFAHopcroft(false).to_map();
  myassert(simp_automaton.size() == 2);
  myassert(simp_automaton.at(0).size() == 8192);
  myassert(simp_automaton.at(32).size() == 32768-8192);
  const auto full_automaton = tft.MinimizeDFAHopcroft(true).to_map();
  myassert(full_automaton.size() == 2);
  myassert(full_automaton.at(0).size() == 8192);
  myassert(full_automaton.at(32).size() == 32768-8192);
}

void test_WSLS() {
  const StrategyN3M5 wsls = StrategyN3M5::WSLS();
  myassert(wsls.IsDefensibleDFA() == false);
  myassert(wsls.IsEfficient() == true);
  myassert(wsls.IsEfficientTopo() == true);
  auto dests = wsls.DestsOfITG();
  for (int i: dests) { myassert(i == 0); } // all goes to cccccc

  auto stat = wsls.StationaryState(0.0001);
  myassert(stat[0] > 0.99);

  myassert(wsls.IsDistinguishable() == true);
  myassert(wsls.IsDistinguishableTopo() == true);

  const auto simp_automaton = wsls.MinimizeDFAHopcroft(false).to_map();
  myassert(simp_automaton.size() == 2);
  myassert(simp_automaton.at(0).size() == 8192);
  myassert(simp_automaton.at(1).size() == 32768-8192);
  const auto full_automaton = wsls.MinimizeDFAHopcroft(true).to_map();
  myassert(full_automaton.size() == 2);
  myassert(full_automaton.at(0).size() == 8192);
  myassert(full_automaton.at(1).size() == 32768-8192);
}


void test_m3_FUSS() {
  const StrategyN3M5 fuss = StrategyN3M5::FUSS_m3();

  const auto simp_automaton = fuss.MinimizeDFAHopcroft(false);
  myassert(simp_automaton.size() == 12);
  // myassert(simp_automaton.to_map() == fuss.MinimizeDFA(false).to_map());

  const auto full_automaton = fuss.MinimizeDFAHopcroft(true);
  myassert(full_automaton.size() == 26);

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
}

void test_AON5() {
  const StrategyN3M5 aon5 = StrategyN3M5::AON(5);
  std::cerr << aon5 << std::endl;

  const auto simp_automaton = aon5.MinimizeDFAHopcroft(false);
  std::cerr << "simplified automaton: " << simp_automaton << std::endl;
  myassert(simp_automaton.size() == 6 );
  const auto full_a = aon5.MinimizeDFAHopcroft(true).to_map();
  myassert(full_a.size() == 6);

  myassert(aon5.IsDefensibleDFA() == false);
  myassert(aon5.IsEfficient() == true);
  myassert(aon5.IsEfficientTopo() == true);

  auto stat = aon5.StationaryState(0.00001);
  myassert(stat[0] > 0.99);
  std::cerr << "stationary state" << std::endl;
  for (size_t i = 0; i < stat.size(); i++) {
    if (stat[i] > 0.05) { std::cerr << StateN3M5(i).ToString() << " : " << stat[i] << std::endl; }
  }

  myassert(aon5.IsDistinguishable() == true);
  myassert(aon5.IsDistinguishableTopo() == true);
}


void test_CAPRI3() {
  const StrategyN3M5 capri = StrategyN3M5::CAPRI3();

  {
    uint64_t i = StateN3M5("cdccc_dcccc_ccccc").ID();
    auto t = capri.TraceStates(i);
    for (uint64_t i: t) {
      std::cerr << StateN3M5(i) << ", "; }
    std::cerr << std::endl;
  }

  auto dests = capri.DestsOfITG();
  for (uint64_t d : dests) {
    if (d != 0 && d != 32767) myassert(false);
  }

  auto stat = capri.StationaryState(0.00001);
  myassert(stat[0] > 0.99);
  std::cerr << "stationary state" << std::endl;
  for (size_t i = 0; i < stat.size(); i++) {
    if (stat[i] > 0.05) { std::cerr << StateN3M5(i).ToString() << " : " << stat[i] << std::endl; }
  }

  myassert(capri.IsEfficient() == true);
  myassert(capri.IsEfficientTopo() == true);

  myassert(capri.IsDefensibleDFA() == true);

  myassert(capri.IsDistinguishable());
  myassert(capri.IsDistinguishableTopo());

  const auto simp_automaton = capri.MinimizeDFAHopcroft(false);
  std::cerr << "simp_automaton: " << simp_automaton.size() << std::endl;
  myassert(simp_automaton.size() == 379);

  const auto full_automaton = capri.MinimizeDFAHopcroft(true);
  std::cerr << "full_automaton: " << full_automaton.size() << std::endl;
  myassert(full_automaton.size() == 1446);
}

void test_sCAPRI3() {
  StrategyN3M5 scapri = StrategyN3M5::sCAPRI3();

  auto dests = scapri.DestsOfITG();
  for (uint64_t d : dests) {
    if (d != 0 && d != 32767) myassert(false);
  }

  auto stat = scapri.StationaryState(0.00001);
  myassert(stat[0] > 0.99);
  std::cerr << "stationary state" << std::endl;
  for (size_t i = 0; i < stat.size(); i++) {
    if (stat[i] > 0.05) { std::cerr << StateN3M5(i).ToString() << " : " << stat[i] << std::endl; }
  }

  myassert(scapri.IsEfficient() == true);
  myassert(scapri.IsEfficientTopo() == true);

  myassert(scapri.IsDistinguishable() == false);
  myassert(scapri.IsDistinguishableTopo() == false);

  myassert(scapri.IsDefensibleDFA() == true);

  const auto simp_automaton = scapri.MinimizeDFAHopcroft(false);
  std::cerr << "simp_automaton: " << simp_automaton.size() << std::endl;
  myassert(simp_automaton.size() == 49);
  const auto full_automaton = scapri.MinimizeDFAHopcroft(false);
  std::cerr << "full_automaton: " << full_automaton.size() << std::endl;
  myassert(full_automaton.size() == 49);
}


int main() {
  std::cout << "Testing StrategyN3M5 class" << std::endl;

  test_State();
  std::cerr << "Testing AllC" << std::endl;
  test_AllC();
  std::cerr << "Testing AllD" << std::endl;
  test_AllD();
  std::cerr << "Testing TFT" << std::endl;
  test_TFT();
  std::cerr << "Testing WSLS" << std::endl;
  test_WSLS();
  std::cerr << "Testing FUSS" << std::endl;
  test_m3_FUSS();
  std::cerr << "Testing AON5" << std::endl;
  test_AON5();
  std::cerr << "Testing CAPRI3" << std::endl;
  test_CAPRI3();
  std::cerr << "Testing sCAPRI3" << std::endl;
  test_sCAPRI3();
  return 0;
}

