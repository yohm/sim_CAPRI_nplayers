#include <iostream>
#include <cassert>
#include "StrategyN2M3.hpp"

#define myassert(x) do {                              \
if (!(x)) {                                           \
  printf("Assertion failed: %s, file %s, line %d\n"   \
         , #x, __FILE__, __LINE__);                   \
  exit(1);                                            \
  }                                                   \
} while (0)

void test_State() {
  StateN2M3 s("ddcdcd");
  myassert(s.a_3 == D);
  myassert(s.a_2 == D);
  myassert(s.a_1 == C);
  myassert(s.b_3 == D);
  myassert(s.b_2 == C);
  myassert(s.b_1 == D);

  uint64_t id = s.ID();
  myassert(id == 53);

  myassert(s == StateN2M3(id));

  myassert(s.NextState(D, C) == StateN2M3("dcdcdc"));

  myassert(s.RelativePayoff() == -1);
  myassert(StateN2M3("ddcddc").RelativePayoff() == 0);
  myassert(StateN2M3("ccdcdc").RelativePayoff() == 1);

  myassert(StateN2M3("cdcdcd").SwapAB() == StateN2M3("dcdcdc"));

  auto noised = StateN2M3("ddccdc").NoisedStates();
  myassert(noised[0] == StateN2M3("dddcdc"));
  myassert(noised[1] == StateN2M3("ddccdd"));

  auto prev = StateN2M3("ddccdc").PossiblePrevStates();
  myassert(prev[0] == StateN2M3("cddccd"));
  myassert(prev[1] == StateN2M3("cdddcd"));
  myassert(prev[2] == StateN2M3("dddccd"));
  myassert(prev[3] == StateN2M3("ddddcd"));
}

void test_AllC() {
  const std::array<Action, 64> acts = {
      C, C, C, C, D, D, D, D,
      C, C, C, C, D, D, D, D,
      C, C, C, C, D, D, D, D,
      C, C, C, C, D, D, D, D,
      C, C, C, C, D, D, D, D,
      C, C, C, C, D, D, D, D,
      C, C, C, C, D, D, D, D,
      C, C, C, C, D, D, D, D
  };
  StrategyN2M3 s1(acts);
  myassert(s1.actions[0] == C);
  myassert(s1.actions[7] == D);
  myassert(s1.actions[59] == C);
  myassert(s1.actions[63] == D);

  std::string bits("ccccddddccccddddccccddddccccddddccccddddccccddddccccddddccccdddd");
  myassert(s1.ToString() == bits);
  myassert(s1 == StrategyN2M3(bits.c_str()));

  myassert(s1.ActionAt(StateN2M3("cccccc")) == C);
  myassert(s1.ActionAt("ddddcc") == D);  // implicit conversion

  myassert(s1.IsDefensible());

  {
    StrategyN2M3 alld("dddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddd");
    myassert(alld.IsDefensible() == true);
    myassert(alld.IsDefensibleDFA() == true);
    myassert(alld.IsEfficient() == false);
    myassert(alld.IsEfficientTopo() == false);
    auto dests = alld.DestsOfITG();
    for (int i: dests) { myassert(i == 63); } // all goes to dddddd

    auto stat = alld.StationaryState(0.001);
    for (int i = 0; i < 63; i++) { myassert(stat[i] < 0.01); }
    myassert(stat[63] > 0.99);

    myassert(alld.IsDistinguishable() == true);
    myassert(alld.IsDistinguishableTopo() == true);

    const auto simp_automaton = alld.MinimizeDFA(false).to_map();
    myassert(simp_automaton.size() == 1);
    const auto full_automaton = alld.MinimizeDFA(true).to_map();
    myassert(full_automaton.size() == 1);

    const auto simp_a = alld.MinimizeDFAHopcroft(false).to_map();
    myassert(simp_a.size() == 1);
    const auto full_a = alld.MinimizeDFAHopcroft(true).to_map();
    myassert(full_a.size() == 1);

  }
  {
    StrategyN2M3 allc("cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc");
    myassert(allc.IsDefensible() == false);
    myassert(allc.IsDefensibleDFA() == false);
    myassert(allc.IsEfficient() == true);
    myassert(allc.IsEfficientTopo() == true);
    auto dests = allc.DestsOfITG();
    for (int i: dests) { myassert(i == 0); } // all goes to cccccc

    auto stat = allc.StationaryState(0.001);
    for (int i = 1; i < 64; i++) { myassert(stat[i] < 0.01); }
    myassert(stat[0] > 0.99);

    myassert(allc.IsDistinguishable() == false);
    myassert(allc.IsDistinguishableTopo() == false);

    const auto simp_automaton = allc.MinimizeDFA(false).to_map();
    myassert(simp_automaton.size() == 1);
    const auto full_automaton = allc.MinimizeDFA(true).to_map();
    myassert(full_automaton.size() == 1);

    const auto simp_a = allc.MinimizeDFAHopcroft(false).to_map();
    myassert(simp_a.size() == 1);
    const auto full_a = allc.MinimizeDFAHopcroft(true).to_map();
    myassert(full_a.size() == 1);
  }
  {
    StrategyN2M3 tft("cdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcd");
    myassert(tft.IsDefensible() == true);
    myassert(tft.IsDefensibleDFA() == true);
    myassert(tft.IsEfficient() == false);
    myassert(tft.IsEfficientTopo() == false);
    auto dests = tft.DestsOfITG();
    for (int i: dests) {
      myassert(i == 0 || i == 63 || StateN2M3("cdcdcd").ID());
    } // all goes to either cccccc, dddddd, cdcdcd

    auto stat = tft.StationaryState(0.001);
    myassert(abs(stat[0] - 0.25) < 0.01);
    myassert(abs(stat[21] - 0.25) < 0.01);
    myassert(abs(stat[42] - 0.25) < 0.01);
    myassert(abs(stat[63] - 0.25) < 0.01);

    myassert(tft.IsDistinguishable() == false);
    myassert(tft.IsDistinguishableTopo() == false);

    const auto simp_automaton = tft.MinimizeDFA(false).to_map();
    myassert(simp_automaton.size() == 2);
    myassert(simp_automaton.at(0).size() == 32);
    myassert(simp_automaton.at(1).size() == 32);
    const auto full_automaton = tft.MinimizeDFA(true).to_map();
    myassert(full_automaton.size() == 2);
    myassert(full_automaton.at(0).size() == 32);
    myassert(full_automaton.at(1).size() == 32);

    const auto simp_a = tft.MinimizeDFAHopcroft(false).to_map();
    myassert(simp_a.size() == 2);
    myassert(simp_automaton == simp_a);
    const auto full_a = tft.MinimizeDFAHopcroft(true).to_map();
    myassert(full_automaton == full_a);
  }
  {
    StrategyN2M3 wsls("cdcdcdcddcdcdcdccdcdcdcddcdcdcdccdcdcdcddcdcdcdccdcdcdcddcdcdcdc");
    myassert(wsls.IsDefensible() == false);
    myassert(wsls.IsDefensibleDFA() == false);
    myassert(wsls.IsEfficient() == true);
    myassert(wsls.IsEfficientTopo() == true);
    auto dests = wsls.DestsOfITG();
    for (int i: dests) { myassert(i == 0); } // all goes to cccccc

    auto stat = wsls.StationaryState(0.001);
    for (int i = 1; i < 64; i++) { myassert(stat[i] < 0.01); }
    myassert(stat[0] > 0.99);

    myassert(wsls.IsDistinguishable() == true);
    myassert(wsls.IsDistinguishableTopo() == true);

    const auto simp_automaton = wsls.MinimizeDFA(false).to_map();
    myassert(simp_automaton.size() == 2);
    myassert(simp_automaton.at(0).size() == 32);
    myassert(simp_automaton.at(1).size() == 32);
    const auto full_automaton = wsls.MinimizeDFA(true).to_map();
    myassert(full_automaton.size() == 2);
    myassert(full_automaton.at(0).size() == 32);
    myassert(full_automaton.at(1).size() == 32);

    const auto simp_a = wsls.MinimizeDFAHopcroft(false).to_map();
    myassert(simp_a.size() == 2);
    const auto full_a = wsls.MinimizeDFAHopcroft(true).to_map();
    myassert(full_automaton == full_a);
  }
  {
    StrategyN2M3 tf2t("cccdcccdcccdcccdcccdcccdcccdcccdcccdcccdcccdcccdcccdcccdcccdcccd"); // tf2t
    myassert(tf2t.IsDefensible() == false);
    myassert(tf2t.IsDefensibleDFA() == false);
    myassert(tf2t.IsEfficient() == true);
    myassert(tf2t.IsEfficientTopo() == true);
    auto dests = tf2t.DestsOfITG();
    for (int i: dests) { myassert(i == 0 || i == 63); }

    auto stat = tf2t.StationaryState(0.001);
    myassert(stat[0] > 0.99);
    for (int i = 1; i < 64; i++) { myassert(stat[i] < 0.01); }

    myassert(tf2t.IsDistinguishable() == false);
    myassert(tf2t.IsDistinguishableTopo() == false);

    const auto simp_automaton = tf2t.MinimizeDFA(false).to_map();
    myassert(simp_automaton.size() == 3);
    myassert(simp_automaton.at(0).size() == 32);
    myassert(simp_automaton.at(1) == std::set<size_t>({1, 5, 9, 13, 17, 21, 25, 29, 33, 37, 41, 45, 49, 53, 57, 61}));
    myassert(simp_automaton.at(3) == std::set<size_t>({3, 7, 11, 15, 19, 23, 27, 31, 35, 39, 43, 47, 51, 55, 59, 63}));

    const auto full_automaton = tf2t.MinimizeDFA(true).to_map();
    myassert(full_automaton.size() == 3);
    myassert(full_automaton.at(0).size() == 32);
    myassert(full_automaton.at(1) == std::set<size_t>({1, 5, 9, 13, 17, 21, 25, 29, 33, 37, 41, 45, 49, 53, 57, 61}));
    myassert(full_automaton.at(3) == std::set<size_t>({3, 7, 11, 15, 19, 23, 27, 31, 35, 39, 43, 47, 51, 55, 59, 63}));

    const auto simp_a = tf2t.MinimizeDFAHopcroft(false).to_map();
    myassert(simp_a.size() == 3);
    myassert(simp_automaton == simp_a);
    const auto full_a = tf2t.MinimizeDFAHopcroft(true).to_map();
    myassert(full_automaton == full_a);
  }

  {
    StrategyN2M3 s("ccddcdddccccdccdcdddddccdccccccdcdccccdcdccddccdcccdddccdccccccd");
    myassert(s.IsEfficient() == true);
    auto stat = s.StationaryState(0.0001);
    myassert(s.IsEfficientTopo() == true);
  }
}

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

  const auto simp_a = tft_atft.MinimizeDFAHopcroft(false).to_map();
  myassert(simp_a.size() == 4);
  myassert(simp_automaton == simp_a);
  const auto full_a = tft_atft.MinimizeDFAHopcroft(true).to_map();
  myassert(full_auto == full_a);
}

void test_AON() {
  for (size_t n = 1; n <= 3; n++) {
    StrategyN2M3 aon = StrategyN2M3::AON(n);
    myassert(aon.IsDefensible() == false);
    myassert(aon.IsDefensibleDFA() == false);

    myassert(aon.IsEfficient());
    myassert(aon.IsEfficientTopo());

    myassert(aon.IsDistinguishable());
    myassert(aon.IsDistinguishableTopo());
  }
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

  std::cerr << capri << std::endl;

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

  const auto simp_a = capri.MinimizeDFAHopcroft(false).to_map();
  myassert(simp_a.size() == 7);
  myassert(simp_auto == simp_a);
  const auto full_a = capri.MinimizeDFAHopcroft(true).to_map();
  myassert(full_auto == full_a);
}


void test_CAPRI2() {
  StrategyN2M3 capri2 = StrategyN2M3::CAPRI2();
  std::cerr << "CAPRI-2" << std::endl;
  std::cerr << capri2 << std::endl;

  std::cerr << "difference of CAPRI-2 from CAPRI" << std::endl;
  StrategyN2M3 capri = StrategyN2M3::CAPRI();
  for (size_t i = 0; i < 64; i++) {
    if (capri.ActionAt(i) != capri2.ActionAt(i)) {
      std::cerr << StateN2M3(i) << " | " << capri2.ActionAt(i) << " : " << capri.ActionAt(i) << std::endl;
    }
  }

  auto stationary = capri2.StationaryState(0.0001, &capri2);
  double coop_level = 0.0;
  for (size_t i = 0; i < 64; i++) {
    StateN2M3 s(i);
    double p = stationary.at(i);
    int n_c = 0;
    if (s.a_1 == C) n_c++;
    if (s.b_1 == C) n_c++;
    coop_level += n_c * 0.5 * p;
  }
  std::cerr << "CAPRI2 coop_level: " << coop_level << std::endl;

  myassert(capri2.IsEfficient());
  myassert(capri2.IsEfficientTopo());
  myassert(capri2.IsDefensible());
  myassert(capri2.IsDefensibleDFA());
  myassert(capri2.IsDistinguishable());
  myassert(capri2.IsDistinguishableTopo());

  const auto simp_auto = capri2.MinimizeDFA(false).to_map();
  myassert(simp_auto.size() == 10);
  const auto full_auto = capri2.MinimizeDFA(true).to_map();
  myassert(full_auto.size() == 18);

  const auto simp_a = capri2.MinimizeDFAHopcroft(false).to_map();
  myassert(simp_auto == simp_a);
  const auto full_a = capri2.MinimizeDFAHopcroft(true).to_map();
  myassert(full_auto == full_a);

  StrategyN2M3 wsls = StrategyN2M3("cdcdcdcddcdcdcdccdcdcdcddcdcdcdccdcdcdcddcdcdcdccdcdcdcddcdcdcdc");
  auto ss = capri2.StationaryState2(0.0001, &wsls);
  std::cerr << ss[0] << std::endl;
  myassert(ss[0] < 0.5); // distinguishable against WSLS
}


void test_sCAPRI2() {
  StrategyN2M3 capri = StrategyN2M3::CAPRI();

  StrategyN2M3 scapri2 = StrategyN2M3::sCAPRI2();
  std::cerr << "scapri2" << std::endl << scapri2 << std::endl;
  std::cerr << "difference of sCAPRI-2 from CAPRI" << std::endl;
  for (size_t i = 0; i < 64; i++) {
    if (capri.ActionAt(i) != scapri2.ActionAt(i)) {
      std::cerr << StateN2M3(i) << " | " << scapri2.ActionAt(i) << " : " << capri.ActionAt(i) << std::endl;
    }
  }

  myassert(scapri2.IsEfficient());
  myassert(scapri2.IsEfficientTopo());
  myassert(scapri2.IsDefensible());
  myassert(scapri2.IsDefensibleDFA());
  myassert(scapri2.IsDistinguishable() == false);
  myassert(scapri2.IsDistinguishableTopo() == false);

  const auto simp_auto = scapri2.MinimizeDFA(false).to_map();
  const auto simp_a = scapri2.MinimizeDFAHopcroft(false).to_map();
  myassert(simp_auto.size() == 7);
  myassert(simp_auto == simp_a);

  const auto full_auto = scapri2.MinimizeDFA(true).to_map();
  myassert(full_auto.size() == 10);
  const auto full_a = scapri2.MinimizeDFAHopcroft(true).to_map();
  myassert(full_auto == full_a);
}

void test_EfficiencyDefensible() {
  StrategyN2M3 s1("cdddddddcdcdddcddcddcdddddccdcddddcdcdddddddddcdddddddcddddcdddd");  // efficient and defensible
  myassert(s1.IsEfficient());
  myassert(s1.IsDefensible());
  myassert(s1.IsDefensibleDFA());
  // auto stat = s1.StationaryState(0.0001);
  // for(int i=0; i<64; i++) { myassert(stat[i] < 0.01); }
  const auto simp_auto = s1.MinimizeDFA(false).to_map();
  const auto simp_a = s1.MinimizeDFAHopcroft(false).to_map();
  myassert(simp_auto == simp_a);
  const auto full_auto = s1.MinimizeDFA(true).to_map();
  const auto full_a = s1.MinimizeDFAHopcroft(true).to_map();
  myassert(full_auto == full_a);
}

void PrintCAPRIs() {
  StrategyN2M3 capri("cdddcdddcdcddddddcddcdddddddddddcdcdcdcdddddddddddddcdccddddddcd");
  const size_t N = 64;
  std::array<Action,N> X,Y;
  StrategyN2M3 capri2 = StrategyN2M3::CAPRI2();
  StrategyN2M3 scapri2 = StrategyN2M3::sCAPRI2();

  for (int i = 0; i < 8; i++) {
    std::cout << i;
    for (int j = 0; j < 8; j++) {
      StateN2M3 s(i*8 + j);
      // std::cout << " & $" << capri.ActionAt(s) << "," << scapri2.ActionAt(s) << "," << capri2.ActionAt(s) << "$";
      std::cout << " & $" << capri2.ActionAt(s) << "$";
    }
    std::cout << " \\\\" << std::endl;
  }
}

int main() {
  std::cout << "Testing StrategyN2M3 class" << std::endl;

  test_State();
  test_AllC();
  test_EfficiencyDefensible();
  test_TFTATFT();
  test_AON();
  test_CAPRI();
  test_CAPRI2();
  test_sCAPRI2();
  PrintCAPRIs();

  return 0;
}

