#include <iostream>
#include <cassert>
#include "StrategyN3M5.hpp"

void test_State() {
  StateN3M5 s("dcdcd_ddccd_ccccc");
  assert(s.ha[0] == D && s.ha[1] == C && s.ha[2] == D && s.ha[3] == C && s.ha[4] == D);
  assert(s.hb[0] == D && s.hb[1] == D && s.hb[2] == C && s.hb[3] == C && s.hb[4] == D);
  assert(s.hc[0] == C && s.hc[1] == C && s.hc[2] == C && s.hc[3] == C && s.hc[4] == C);

  assert(s.ToString() == "dcdcd_ddccd_ccccc");

  assert(s.NextState(D, C, D) == StateN3M5("ddcdc_cddcc_dcccc"));

  assert(s.RelativePayoff(true) == 0);
  assert(s.RelativePayoff(false) == 1);

  assert(s.ID() == 0b000001001110101);

  assert(StateN3M5(0b111110101011100) == StateN3M5("ccddd_cdcdc_ddddd"));

  assert(StateN3M5("ddddd_dcdcd_ccccc").StateFromB() == StateN3M5("dcdcd_ccccc_ddddd"));
  assert(StateN3M5("ddddd_dcdcd_ccccc").StateFromC() == StateN3M5("ccccc_ddddd_dcdcd"));

  auto noised = StateN3M5("ddddd_dcdcd_ccccc").NoisedStates();
  assert(noised.size() == 3);
  assert(std::find(noised.cbegin(), noised.cend(), StateN3M5("cdddd_dcdcd_ccccc")) != noised.end());
  assert(std::find(noised.cbegin(), noised.cend(), StateN3M5("ddddd_ccdcd_ccccc")) != noised.end());
  assert(std::find(noised.cbegin(), noised.cend(), StateN3M5("ddddd_dcdcd_dcccc")) != noised.end());

  auto prev = StateN3M5("ddddd_dcdcd_ccccc").PossiblePrevStates();
  assert(prev.size() == 8);
  std::set<std::string> prev_set;
  for (const auto &st: prev) { prev_set.insert(st.ToString()); }
  std::set<std::string> expected = {"ddddc_cdcdc_ccccc", "ddddc_cdcdc_ccccd", "ddddc_cdcdd_ccccc", "ddddc_cdcdd_ccccd",
                                    "ddddd_cdcdc_ccccc", "ddddd_cdcdc_ccccd", "ddddd_cdcdd_ccccc", "ddddd_cdcdd_ccccd"};
  assert(prev_set == expected);

  auto s1 = StateN3M5("ddddd_dcdcd_ccccc");
  assert(s1.NumDiffInT1(StateN3M5("ddddd_dcdcd_ccccc")) == 0);
  assert(s1.NumDiffInT1(StateN3M5("cdddd_dcdcd_ccccc")) == 1);
  assert(s1.NumDiffInT1(StateN3M5("cdddd_ccdcd_ccccc")) == 2);
  assert(s1.NumDiffInT1(StateN3M5("cdddd_ccdcd_dcccc")) == 3);
  assert(s1.NumDiffInT1(StateN3M5("ddddd_dcdcd_ccccd")) == -1);
}

void test_Strategy() {

  {
    std::bitset<32768> allc_b(0ull);
    StrategyN3M5 allc(allc_b);

    std::string s = allc.ToString();
    assert(s.size() == 32768);
    bool b = true;
    for (char c: s) { if(c != '0') { b = false; break; } }
    assert(b);

    assert(allc.ActionAt(StateN3M5(0)) == C);
    assert(allc.ActionAt(StateN3M5(999)) == C);
    assert(allc.ActionAt(StateN3M5(9999)) == C);

    assert(allc.IsEfficientTopo());
    // assert(allc.IsEfficient());

    // UnionFind uf = allc.MinimizeDFA(false);
    // assert(uf.to_map().size() == 1);
  }

  /*
  {
    StrategyN2M3 alld("dddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddd");
    assert(alld.IsDefensible() == true);
    assert(alld.IsDefensibleDFA() == true);
    assert(alld.IsEfficient() == false);
    assert(alld.IsEfficientTopo() == false);
    auto dests = alld.DestsOfITG();
    for (int i: dests) { assert(i == 63); } // all goes to dddddd

    auto stat = alld.StationaryState(0.001);
    for (int i = 0; i < 63; i++) { assert(stat[i] < 0.01); }
    assert(stat[63] > 0.99);

    assert(alld.IsDistinguishable() == true);
    assert(alld.IsDistinguishableTopo() == true);

    const auto simp_automaton = alld.MinimizeDFA(false).to_map();
    assert(simp_automaton.size() == 1);
    const auto full_automaton = alld.MinimizeDFA(true).to_map();
    assert(full_automaton.size() == 1);
  }
  {
    StrategyN2M3 allc("cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc");
    assert(allc.IsDefensible() == false);
    assert(allc.IsDefensibleDFA() == false);
    assert(allc.IsEfficient() == true);
    assert(allc.IsEfficientTopo() == true);
    auto dests = allc.DestsOfITG();
    for (int i: dests) { assert(i == 0); } // all goes to cccccc

    auto stat = allc.StationaryState(0.001);
    for (int i = 1; i < 64; i++) { assert(stat[i] < 0.01); }
    assert(stat[0] > 0.99);

    assert(allc.IsDistinguishable() == false);
    assert(allc.IsDistinguishableTopo() == false);

    const auto simp_automaton = allc.MinimizeDFA(false).to_map();
    assert(simp_automaton.size() == 1);
    const auto full_automaton = allc.MinimizeDFA(true).to_map();
    assert(full_automaton.size() == 1);
  }
  {
    StrategyN2M3 tft("cdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcd");
    assert(tft.IsDefensible() == true);
    assert(tft.IsDefensibleDFA() == true);
    assert(tft.IsEfficient() == false);
    assert(tft.IsEfficientTopo() == false);
    auto dests = tft.DestsOfITG();
    for (int i: dests) {
      assert(i == 0 || i == 63 || StateN2M3("cdcdcd").ID());
    } // all goes to either cccccc, dddddd, cdcdcd

    auto stat = tft.StationaryState(0.001);
    assert(abs(stat[0] - 0.25) < 0.01);
    assert(abs(stat[21] - 0.25) < 0.01);
    assert(abs(stat[42] - 0.25) < 0.01);
    assert(abs(stat[63] - 0.25) < 0.01);

    assert(tft.IsDistinguishable() == false);
    assert(tft.IsDistinguishableTopo() == false);

    const auto simp_automaton = tft.MinimizeDFA(false).to_map();
    assert(simp_automaton.size() == 2);
    assert(simp_automaton.at(0).size() == 32);
    assert(simp_automaton.at(1).size() == 32);
    const auto full_automaton = tft.MinimizeDFA(true).to_map();
    assert(full_automaton.size() == 2);
    assert(full_automaton.at(0).size() == 32);
    assert(full_automaton.at(1).size() == 32);
  }
  {
    StrategyN2M3 wsls("cdcdcdcddcdcdcdccdcdcdcddcdcdcdccdcdcdcddcdcdcdccdcdcdcddcdcdcdc");
    assert(wsls.IsDefensible() == false);
    assert(wsls.IsDefensibleDFA() == false);
    assert(wsls.IsEfficient() == true);
    assert(wsls.IsEfficientTopo() == true);
    auto dests = wsls.DestsOfITG();
    for (int i: dests) { assert(i == 0); } // all goes to cccccc

    auto stat = wsls.StationaryState(0.001);
    for (int i = 1; i < 64; i++) { assert(stat[i] < 0.01); }
    assert(stat[0] > 0.99);

    assert(wsls.IsDistinguishable() == true);
    assert(wsls.IsDistinguishableTopo() == true);

    const auto simp_automaton = wsls.MinimizeDFA(false).to_map();
    assert(simp_automaton.size() == 2);
    assert(simp_automaton.at(0).size() == 32);
    assert(simp_automaton.at(1).size() == 32);
    const auto full_automaton = wsls.MinimizeDFA(true).to_map();
    assert(full_automaton.size() == 2);
    assert(full_automaton.at(0).size() == 32);
    assert(full_automaton.at(1).size() == 32);
  }
  {
    StrategyN2M3 tf2t("cccdcccdcccdcccdcccdcccdcccdcccdcccdcccdcccdcccdcccdcccdcccdcccd"); // tf2t
    assert(tf2t.IsDefensible() == false);
    assert(tf2t.IsDefensibleDFA() == false);
    assert(tf2t.IsEfficient() == true);
    assert(tf2t.IsEfficientTopo() == true);
    auto dests = tf2t.DestsOfITG();
    for (int i: dests) { assert(i == 0 || i == 63); }

    auto stat = tf2t.StationaryState(0.001);
    assert(stat[0] > 0.99);
    for (int i = 1; i < 64; i++) { assert(stat[i] < 0.01); }

    assert(tf2t.IsDistinguishable() == false);
    assert(tf2t.IsDistinguishableTopo() == false);

    const auto simp_automaton = tf2t.MinimizeDFA(false).to_map();
    assert(simp_automaton.size() == 3);
    assert(simp_automaton.at(0).size() == 32);
    assert(simp_automaton.at(1) == std::set<size_t>({1, 5, 9, 13, 17, 21, 25, 29, 33, 37, 41, 45, 49, 53, 57, 61}));
    assert(simp_automaton.at(3) == std::set<size_t>({3, 7, 11, 15, 19, 23, 27, 31, 35, 39, 43, 47, 51, 55, 59, 63}));

    const auto full_automaton = tf2t.MinimizeDFA(true).to_map();
    assert(full_automaton.size() == 3);
    assert(full_automaton.at(0).size() == 32);
    assert(full_automaton.at(1) == std::set<size_t>({1, 5, 9, 13, 17, 21, 25, 29, 33, 37, 41, 45, 49, 53, 57, 61}));
    assert(full_automaton.at(3) == std::set<size_t>({3, 7, 11, 15, 19, 23, 27, 31, 35, 39, 43, 47, 51, 55, 59, 63}));
  }

  {
    StrategyN2M3 s("ccddcdddccccdccdcdddddccdccccccdcdccccdcdccddccdcccdddccdccccccd");
    assert(s.IsEfficient() == true);
    auto stat = s.StationaryState(0.0001);
    assert(s.IsEfficientTopo() == true);
  }
   */
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

  assert(tft_atft.ToString() == "cdcdcdcddccddccdcdcccdccdccddccdcdcdcdcddccddccdcdcccdccdccddccd");

  assert(tft_atft.IsDefensible());
  assert(tft_atft.IsDefensibleDFA());
  assert(tft_atft.IsEfficient());

  assert(tft_atft.IsDistinguishable());
  assert(tft_atft.IsDistinguishableTopo());

  const auto simp_automaton = tft_atft.MinimizeDFA(false).to_map();
  assert(simp_automaton.size() == 4);
  assert(simp_automaton.at(0).size() == 28);  // TFT-c
  assert(simp_automaton.at(1).size() == 20);  // TFT-d
  assert(simp_automaton.at(8) == std::set<size_t>({8, 12, 40, 44, 24, 28, 56, 60}));  // ATFT-d
  assert(simp_automaton.at(9) == std::set<size_t>({9, 13, 41, 45, 25, 29, 57, 61}));  // ATFT-c

  const auto full_auto = tft_atft.MinimizeDFA(true).to_map();
  assert(full_auto.size() == 6);
  assert(full_auto.at(0).size() == 24);  // TFT-c
  assert(full_auto.at(1).size() == 12);  // TFT-d
  assert(full_auto.at(8) == std::set<size_t>({8, 12, 40, 44, 24, 28, 56, 60}));  // ATFT-d
  assert(full_auto.at(9) == std::set<size_t>({9, 13, 41, 45, 25, 29, 57, 61}));  // ATFT-c
  assert(full_auto.at(19) == std::set<size_t>({19, 23, 51, 55}));  // TFT-c-2
  assert(full_auto.at(11) == std::set<size_t>({11, 15, 43, 47, 27, 31, 59, 63}));  // TFT-d-2
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

  assert(capri.IsDefensible());
  assert(capri.IsDefensibleDFA());
  assert(capri.IsEfficient());
  assert(capri.IsEfficientTopo());

  assert(capri.IsDistinguishable());
  assert(capri.IsDistinguishableTopo());

  { // distinguishable against WSLS
    StrategyN2M3 wsls("cdcdcdcddcdcdcdccdcdcdcddcdcdcdccdcdcdcddcdcdcdccdcdcdcddcdcdcdc");
    const auto stat = capri.StationaryState(1.0e-6, &wsls);
    assert(std::abs(stat[61] - 0.5) < 0.01);  // ddd,dcd ~ 0.5
    assert(std::abs(stat[58] - 0.5) < 0.01);  // ddd,cdc ~ 0.5
    assert(stat[0] < 0.01);  // ccc,ccc ~ 0.5
  }

  const auto simp_auto = capri.MinimizeDFA(false).to_map();
  assert(simp_auto.size() == 7);
  const auto full_auto = capri.MinimizeDFA(true).to_map();
  assert(full_auto.size() == 14);
}

void test_EfficiencyDefensible() {
  StrategyN2M3 s1("cdddddddcdcdddcddcddcdddddccdcddddcdcdddddddddcdddddddcddddcdddd");  // efficient and defensible
  assert(s1.IsEfficient());
  assert(s1.IsDefensible());
  assert(s1.IsDefensibleDFA());
  // auto stat = s1.StationaryState(0.0001);
  // for(int i=0; i<64; i++) { assert(stat[i] < 0.01); }
}
 */

int main() {
  std::cout << "Testing StrategyN3M5 class" << std::endl;

  test_State();
  test_Strategy();
  // test_EfficiencyDefensible();
  // test_TFTATFT();
  // test_CAPRI();
  return 0;
}

