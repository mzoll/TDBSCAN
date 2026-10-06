//
// Created by mzoll on 31/08/2026.
//


#include "external/common_clib/InterruptableQueue.hpp"


int main() {
  auto a = new common_clib::threading::InterruptableQueue<int>();

  for (int i=0; i <4 ; i++)
    a->push(i);

  while (! a->empty())
    const auto _ = a->pop();

  delete a;

  return 0;
}
