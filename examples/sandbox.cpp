//
// Created by netsu on 31/08/2026.
//

#include <set>
#include <iostream>

#include "tdbscan_algo/common_defs.h"

// using namespace tdbscan_algo;
using namespace std;

struct inverse_int {
  int value_;
  bool operator<(const inverse_int rhs) const {
    return value_ > rhs.value_;
  }
  inverse_int(int a) : value_(a) {};
};



int main() {



  std::set<ScalarTime_t> s;
  s.insert(ScalarTime_t(1));
  s.insert(ScalarTime_t(9));
  s.insert(ScalarTime_t(5));

  cout << "abc "<< endl;
  for (const auto& a : s)
    cout << a << endl;


  return 0;
}