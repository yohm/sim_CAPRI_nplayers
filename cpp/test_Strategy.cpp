#include <iostream>
#include <cassert>
#include "Strategy.hpp"

void test_State() {
  State s("ddcdcd");
  assert(s.a_3 == D);
  assert(s.a_2 == D);
  assert(s.a_1 == C);
  assert(s.b_3 == D);
  assert(s.b_2 == C);
  assert(s.b_1 == D);

  uint64_t id = s.ID();
  assert( id == 53 );

  assert( s == State(id) );

  assert( s.NextState(D,C) == State("dcdcdc") );

  assert( s.RelativePayoff() == -1 );
  assert( State("ddcddc").RelativePayoff() == 0 );
  assert( State("ccdcdc").RelativePayoff() == 1 );

  assert( State("cdcdcd").SwapAB() == State("dcdcdc") );

  auto noised = State("ddccdc").NoisedStates();
  assert( noised[0] == State("dddcdc") );
  assert( noised[1] == State("ddccdd") );

  auto prev = State("ddccdc").PossiblePrevStates();
  assert( prev[0] == State("cddccd") );
  assert( prev[1] == State("cdddcd") );
  assert( prev[2] == State("dddccd") );
  assert( prev[3] == State("ddddcd") );
}

void test_Strategy() {
  const std::array<Action,64> acts = {
      C,C,C,C,D,D,D,D,
      C,C,C,C,D,D,D,D,
      C,C,C,C,D,D,D,D,
      C,C,C,C,D,D,D,D,
      C,C,C,C,D,D,D,D,
      C,C,C,C,D,D,D,D,
      C,C,C,C,D,D,D,D,
      C,C,C,C,D,D,D,D
  };
  Strategy s1(acts);
  assert( s1.actions[0] == C );
  assert( s1.actions[7] == D );
  assert( s1.actions[59] == C );
  assert( s1.actions[63] == D );

  std::string bits("ccccddddccccddddccccddddccccddddccccddddccccddddccccddddccccdddd");
  assert( s1.ToString() == bits );
  assert( s1 == Strategy(bits.c_str()) );

  assert( s1.ActionAt(State("cccccc")) == C );
  assert( s1.ActionAt("ddddcc") == D );  // implicit conversion

  assert( s1.IsDefensible() );

  {
    Strategy alld("dddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddd");
    assert( alld.IsDefensible() == true );
    assert( alld.IsEfficient() == false );
    assert( alld.IsEfficientTopo() == false );
    auto dests = alld.DestsOfITG();
    for(int i: dests) { assert( i == 63 ); } // all goes to dddddd

    auto stat = alld.StationaryState(0.001);
    for(int i=0; i<63; i++) { assert(stat[i] < 0.01); }
    assert(stat[63] > 0.99);

    assert( alld.IsDistinguishable() == true );
    assert( alld.IsDistinguishableTopo() == true );

    const auto simp_automaton = alld.MinimizeDFA(false).to_map();
    assert(simp_automaton.size() == 1);
    const auto full_automaton = alld.MinimizeDFA(true).to_map();
    assert(full_automaton.size() == 1);
  }
  {
    Strategy allc("cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc");
    assert( allc.IsDefensible() == false );
    assert( allc.IsEfficient() == true );
    assert( allc.IsEfficientTopo() == true );
    auto dests = allc.DestsOfITG();
    for(int i: dests) { assert( i == 0 ); } // all goes to cccccc

    auto stat = allc.StationaryState(0.001);
    for(int i=1; i<64; i++) { assert(stat[i] < 0.01); }
    assert(stat[0] > 0.99);

    assert( allc.IsDistinguishable() == false );
    assert( allc.IsDistinguishableTopo() == false );

    const auto simp_automaton = allc.MinimizeDFA(false).to_map();
    assert(simp_automaton.size() == 1);
    const auto full_automaton = allc.MinimizeDFA(true).to_map();
    assert(full_automaton.size() == 1);
  }
  {
    Strategy tft("cdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcd");
    assert( tft.IsDefensible() == true );
    assert( tft.IsEfficient() == false );
    assert( tft.IsEfficientTopo() == false );
    auto dests = tft.DestsOfITG();
    for(int i: dests) { assert( i == 0 || i == 63 || State("cdcdcd").ID() ); } // all goes to either cccccc, dddddd, cdcdcd

    auto stat = tft.StationaryState(0.001);
    assert( abs(stat[0]-0.25) < 0.01 );
    assert( abs(stat[21]-0.25) < 0.01 );
    assert( abs(stat[42]-0.25) < 0.01 );
    assert( abs(stat[63]-0.25) < 0.01 );

    assert( tft.IsDistinguishable() == false );
    assert( tft.IsDistinguishableTopo() == false );

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
    Strategy wsls("cdcdcdcddcdcdcdccdcdcdcddcdcdcdccdcdcdcddcdcdcdccdcdcdcddcdcdcdc");
    assert( wsls.IsDefensible() == false );
    assert( wsls.IsEfficient() == true );
    assert( wsls.IsEfficientTopo() == true );
    auto dests = wsls.DestsOfITG();
    for(int i: dests) { assert( i == 0 ); } // all goes to cccccc

    auto stat = wsls.StationaryState(0.001);
    for(int i=1; i<64; i++) { assert(stat[i] < 0.01); }
    assert(stat[0] > 0.99);

    assert( wsls.IsDistinguishable() == true );
    assert( wsls.IsDistinguishableTopo() == true );

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
    Strategy tf2t("cccdcccdcccdcccdcccdcccdcccdcccdcccdcccdcccdcccdcccdcccdcccdcccd"); // tf2t
    assert( tf2t.IsDefensible() == false );
    assert( tf2t.IsEfficient() == true );
    assert( tf2t.IsEfficientTopo() == true );
    auto dests = tf2t.DestsOfITG();
    for(int i: dests) { assert( i == 0 || i == 63 ); }

    auto stat = tf2t.StationaryState(0.001);
    assert(stat[0] > 0.99);
    for(int i=1; i<64; i++) { assert(stat[i] < 0.01); }

    assert( tf2t.IsDistinguishable() == false );
    assert( tf2t.IsDistinguishableTopo() == false );

    const auto simp_automaton = tf2t.MinimizeDFA(false).to_map();
    assert(simp_automaton.size() == 3);
    assert(simp_automaton.at(0).size() == 32);
    assert(simp_automaton.at(1) == std::set<size_t>({1,5,9,13,17,21,25,29,33,37,41,45,49,53,57,61}));
    assert(simp_automaton.at(3) == std::set<size_t>({3,7,11,15,19,23,27,31,35,39,43,47,51,55,59,63}));

    const auto full_automaton = tf2t.MinimizeDFA(true).to_map();
    assert(full_automaton.size() == 3);
    assert(full_automaton.at(0).size() == 32);
    assert(full_automaton.at(1) == std::set<size_t>({1,5,9,13,17,21,25,29,33,37,41,45,49,53,57,61}));
    assert(full_automaton.at(3) == std::set<size_t>({3,7,11,15,19,23,27,31,35,39,43,47,51,55,59,63}));
  }

  {
    Strategy s("ccddcdddccccdccdcdddddccdccccccdcdccccdcdccddccdcccdddccdccccccd");
    assert( s.IsEfficient() == true );
    auto stat = s.StationaryState(0.0001);
    assert( s.IsEfficientTopo() == true );
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
  std::map<int,Action> m = {
      {0, C}, {1, D}, {2, C}, {3, D},
      {8, D}, {9, C}, {10, C}, {11, D},
      {16, C}, {17, D}, {18, C}, {19, C},
      {24, D}, {25, C}, {26, C}, {27, D}
  };
  Strategy tft_atft("cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc");
  int mask = 27; // 011011
  for(int i=0; i<64; i++) {
    int masked = i & mask;
    tft_atft.SetAction(i, m.at(masked) );
  }

  assert( tft_atft.ToString() == "cdcdcdcddccddccdcdcccdccdccddccdcdcdcdcddccddccdcdcccdccdccddccd" );

  assert( tft_atft.IsDefensible() );
  assert( tft_atft.IsEfficient() );

  assert( tft_atft.IsDistinguishable() );
  assert( tft_atft.IsDistinguishableTopo() );

  const auto simp_automaton = tft_atft.MinimizeDFA(false).to_map();
  assert(simp_automaton.size() == 4);
  assert(simp_automaton.at(0).size() == 28);  // TFT-c
  assert(simp_automaton.at(1).size() == 20);  // TFT-d
  assert(simp_automaton.at(8) == std::set<size_t>({8,12,40,44, 24,28,56,60}));  // ATFT-d
  assert(simp_automaton.at(9) == std::set<size_t>({9,13,41,45, 25,29,57,61}));  // ATFT-c

  const auto full_auto = tft_atft.MinimizeDFA(true).to_map();
  assert(full_auto.size() == 6);
  assert(full_auto.at(0).size() == 24);  // TFT-c
  assert(full_auto.at(1).size() == 12);  // TFT-d
  assert(full_auto.at(8) == std::set<size_t>({8,12,40,44, 24,28,56,60}));  // ATFT-d
  assert(full_auto.at(9) == std::set<size_t>({9,13,41,45, 25,29,57,61}));  // ATFT-c
  assert(full_auto.at(19) == std::set<size_t>({19,23,51,55}));  // TFT-c-2
  assert(full_auto.at(11) == std::set<size_t>({11,15,43,47, 27,31,59,63}));  // TFT-d-2
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
  Strategy capri("cdddcdddcdcddddddcddcdddddddddddcdcdcdcdddddddddddddcdccddddddcd");

  assert( capri.IsDefensible() );
  assert( capri.IsEfficient() );
  assert( capri.IsEfficientTopo() );

  assert( capri.IsDistinguishable() );
  assert( capri.IsDistinguishableTopo() );

  { // distinguishable against WSLS
    Strategy wsls("cdcdcdcddcdcdcdccdcdcdcddcdcdcdccdcdcdcddcdcdcdccdcdcdcddcdcdcdc");
    const auto stat = capri.StationaryState(1.0e-6, &wsls);
    assert( std::abs(stat[61]-0.5) < 0.01 );  // ddd,dcd ~ 0.5
    assert( std::abs(stat[58]-0.5) < 0.01 );  // ddd,cdc ~ 0.5
    assert( stat[0] < 0.01 );  // ccc,ccc ~ 0.5
  }

  const auto simp_auto = capri.MinimizeDFA(false).to_map();
  assert(simp_auto.size() == 7);
  const auto full_auto = capri.MinimizeDFA(true).to_map();
  assert(full_auto.size() == 14);
}

void test_EfficiencyDefensible() {
  Strategy s1("cdddddddcdcdddcddcddcdddddccdcddddcdcdddddddddcdddddddcddddcdddd");  // efficient and defensible
  assert( s1.IsEfficient() );
  assert( s1.IsDefensible() );
  // auto stat = s1.StationaryState(0.0001);
  // for(int i=0; i<64; i++) { assert(stat[i] < 0.01); }
}

int main() {
  std::cout << "Testing Strategy class" << std::endl;

  test_State();
  test_Strategy();
  test_EfficiencyDefensible();
  test_TFTATFT();
  test_CAPRI();
  return 0;
}

