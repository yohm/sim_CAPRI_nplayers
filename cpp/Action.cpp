#include <iostream>
#include <cassert>
#include "Action.hpp"

char A2C( Action act ) {
  if(act == C) { return 'c'; }
  else if(act == D) { return 'd'; }
  else if(act == U) { return '_'; }
  else if(act == W) { return '*'; }
  return ' ';
}

Action C2A( char c ) {
  if(c == 'c') { return C; }
  else if(c == 'd') { return D; }
  else if(c == '_') { return U; }
  else if(c == '*') { return W; }
  else {
    assert(false);
    return C;
  }
}

std::ostream &operator<<(std::ostream &os, const Action &act) {
  os << A2C(act);
  return os;
}

