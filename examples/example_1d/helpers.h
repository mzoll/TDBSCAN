//
// Created by netsu on 27/09/2026.
//

#ifndef TDBSCAN__EXAMPLE1D__HELPERS_H
#define TDBSCAN__EXAMPLE1D__HELPERS_H

#include <set>
#include "blib.h"


namespace ex1d {
/**
 * Helper function: calculate the purity, Signal over Noise ratio, of this Blib sample
 * @param blibs the blibs to quantify
 * @return the purity as '#signal / #all'.
 */
double calculate_signal_purity(const std::set<SBlibWithTrace>& blibs) {
  int _signal_count = 0;
  int _noise_count = 0;
  for (const auto& b: blibs) {
    switch (b.origin_) {
      case SBlibWithTrace::SIGNAL:
        _signal_count++;
        break;
      case SBlibWithTrace::NOISE:
        _noise_count++;
        break;
      case SBlibWithTrace::UNKNOWN:
        break;
    }
  }
  return static_cast<double>(_signal_count)/blibs.size();
}

} //namespace ex1d

#endif //TDBSCAN__EXAMPLE1D__HELPERS_H
