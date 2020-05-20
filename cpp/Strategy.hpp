#include <string>
#include <array>
#include <vector>
#include <set>
#include <map>
#include <queue>
#include <stack>
#include <sstream>
#include <cstdint>
#include <ostream>
#include <algorithm>
#include <Eigen/Dense>
#include "Action.hpp"
#include "DirectedGraph.hpp"

#ifndef STRATEGY_HPP
#define STRATEGY_HPP

class Strategy;

class State {
public:
  State(Action _a_3, Action _a_2, Action _a_1, Action _b_3, Action _b_2, Action _b_1):
      a_3(_a_3), a_2(_a_2), a_1(_a_1), b_3(_b_3), b_2(_b_2), b_1(_b_1) {assert(AllCorD());};
  State(uint64_t id):  // upper bit [a_3,a_2,a_1,b_3,b_2,b_1] lower bit
      a_3( ((id>>5)&1)?D:C ), a_2( ((id>>4)&1)?D:C ), a_1( ((id>>3)&1)?D:C ),
      b_3( ((id>>2)&1)?D:C ), b_2( ((id>>1)&1)?D:C ), b_1( ((id>>0)&1)?D:C ) { assert(AllCorD());};
  State(const char str[6]):
      a_3(C2A(str[0])), a_2(C2A(str[1])), a_1(C2A(str[2])),
      b_3(C2A(str[3])), b_2(C2A(str[4])), b_1(C2A(str[5])) {assert(AllCorD());};
  const Action a_3, a_2, a_1, b_3, b_2, b_1;
  bool AllCorD() const {
    return (
        (a_3==D || a_3==C) && (a_2==D || a_2==C) && (a_1==D || a_1==C) &&
        (b_3==D || b_3==C) && (b_2==D || b_2==C) && (b_1==D || b_1==C)
    );
  }

  bool operator==(const State & rhs) const {
    return (a_3==rhs.a_3 && a_2==rhs.a_2 && a_1==rhs.a_1 && b_3==rhs.b_3 && b_2==rhs.b_2 && b_1==rhs.b_1);
  }
  friend std::ostream &operator<<(std::ostream &os, const State &s) {
    os << s.a_3 << s.a_2 << s.a_1 << s.b_3 << s.b_2 << s.b_1;
    return os;
  };

  State NextState(Action act_a, Action act_b) const {
    assert(act_a == C || act_a == D);
    assert(act_b == C || act_b == D);
    return State(a_2,a_1,act_a,b_2,b_1,act_b);
  };

  std::array<State,4> PossiblePrevStates() const {
    std::array<State,4> ans = {
        State(C, a_3, a_2, C, b_3, b_2),
        State(C, a_3, a_2, D, b_3, b_2),
        State(D, a_3, a_2, C, b_3, b_2),
        State(D, a_3, a_2, D, b_3, b_2)
    };
    return ans;
  }

  int RelativePayoff() const {
    if( a_1 == C && b_1 == D ) { return -1; }
    else if( a_1 == D && b_1 == C ) { return  1; }
    else if( a_1 == b_1 ) { return 0; }
    else { assert(false); return -10000; }
  }

  State SwapAB() const { return State(b_3, b_2, b_1, a_3, a_2, a_1); } // state from B's viewpoint

  std::array<State,2> NoisedStates() const {
    Action a_1_n = (a_1==C) ? D : C;
    Action b_1_n = (b_1==C) ? D : C;
    std::array<State,2> ans = { State(a_3,a_2,a_1_n,b_3,b_2,b_1), State(a_3,a_2,a_1,b_3,b_2,b_1_n) };
    return std::move(ans);
  }
  int NumDiffInT1(const State& other) const {
    if(a_3 != other.a_3 || a_2 != other.a_2 || b_3 != other.b_3 || b_2 != other.b_2 ) {
      return -1;
    }
    else {
      int cnt = 0;
      if(a_1 != other.a_1) cnt++;
      if(b_1 != other.b_1) cnt++;
      return cnt;
    }
  }

  uint64_t ID() const {  // ID must be 0~63 integer. AllC: 0, AllD: 63
    uint64_t id = 0;
    if( a_3 == D ) { id += 1 << 5; }
    if( a_2 == D ) { id += 1 << 4; }
    if( a_1 == D ) { id += 1 << 3; }
    if( b_3 == D ) { id += 1 << 2; }
    if( b_2 == D ) { id += 1 << 1; }
    if( b_1 == D ) { id += 1 << 0; }
    return id;
  }
  bool operator<( const State & rhs ) const{
    return ( ID() < rhs.ID() );
  }

};


class Strategy {
public:
  Strategy( const std::array<Action,64>& acts ); // construct a strategy from a list of actions
  Strategy( const char acts[64] );
  std::array<Action,64> actions;

  std::string ToString() const;
  friend std::ostream &operator<<(std::ostream &os, const Strategy &strategy);
  bool operator==(const Strategy & rhs) const {
    for(size_t i=0; i<64; i++) { if(actions[i] != rhs.actions[i]) return false; }
    return true;
  }

  Action ActionAt( const State& s ) const { return actions[s.ID()]; }
  void SetAction( const State& s, Action a ) { assert(actions[s.ID()]==U||actions[s.ID()]==W); actions[s.ID()] = a; d_matrix_ready = false; }
  int NumFixed() const { int c=0; for(auto a: actions) { if(a==C||a==D) c++;} return c; }
  int NumU() const { int c=0; for(auto a: actions) { if(a==U) c++;} return c; }
  uint64_t Size() const { return (1ULL << (64 - NumFixed())); }
  bool IsDefensible();  // check defensibility. If defensible, m_d is also calculated
  bool IsDefensible2() const;  // check defensibility without caching d matrix
  std::array<double,64> StationaryState(double e=0.0001, const Strategy* coplayer = NULL) const; // all actions must be fixed to calculated stationary state
  bool IsEfficient(double e=0.00001, double th=0.95) const {
    return (StationaryState(e)[0]>th);
  } // check efficiency. all actions must be fixed
  bool IsEfficientTopo() const; // check efficiency using ITG
  bool IsDistinguishable(double e=0.00001, double th=0.95) const {
    const Strategy allc("cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc");
    return (StationaryState(e,&allc)[0]<th); };  // check distinguishability against AllC
  bool IsDistinguishableTopo() const; // check distinguishability using the transition graph
  bool SetActionAndRecalcD(const State& s, Action a); // set action[s]=a, and recalculate `m_d`. If not defensible, return false.
  std::vector<State> DanglingStates() const; // states that has incoming links but has no outgoing links (i.e. its action is U)
  std::vector<State> NegativeDanglingStates() const; // states that has incoming links but has no outgoing links (i.e. its action is U). IsDefensible must be called beforehand
  DirectedGraph ITG(bool make_link_UW = true) const;  // construct ITG. When a state is U or W, outgoing links for both {c,d} actions are added when make_link_UW is true.
  std::array<int,64> DestsOfITG() const; // Trace the intra-transition-graph from node i. Destination is stored in i'th element. Undetermined destination is -1
  int NextITGState(const State& s) const; // Trace the intra-transition graph by one step
private:
  typedef std::array<std::array<int8_t,64>,64> d_matrix_t;
  d_matrix_t m_d;
  bool d_matrix_ready;
  std::vector<State> NextPossibleStates(State current) const;
};

#endif //STRATEGY_HPP

