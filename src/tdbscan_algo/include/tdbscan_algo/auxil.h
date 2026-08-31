//
// Created by netsu on 31/08/2026.
//

#ifndef TDBSCAN_AUX_H
#define TDBSCAN_AUX_H

#include <exception>
#include <stdexcept>

/**
 * defines a pure float number on the interval [0. ... 1.]
 * */
class ZeroOne {
  double value_;
public:
  ZeroOne(const double value) : value_(value) {
    if (value < 0. || value > 1.)
      throw std::invalid_argument("value must be within the range [0. ... 1.].");
  }

  operator double() const {return value_;}
  ZeroOne& operator=(const double value) {value_=value; return *this;}
};


#endif //TDBSCAN_AUX_H
