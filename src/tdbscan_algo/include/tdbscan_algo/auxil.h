//
// Created by netsu on 31/08/2026.
//

#ifndef TDBSCAN_AUX_H
#define TDBSCAN_AUX_H

#include <exception>
#include <stdexcept>
#include <random>
#include <sstream>


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



namespace uuid {
static std::random_device              rd;
static std::mt19937                    gen(rd());
static std::uniform_int_distribution<> dis(0, 15);
static std::uniform_int_distribution<> dis2(8, 11);

/// a uuid in string representation
std::string generate_uuid_v4() {
  std::stringstream ss;
  int i;
  ss << std::hex;
  for (i = 0; i < 8; i++) {
    ss << dis(gen);
  }
  ss << "-";
  for (i = 0; i < 4; i++) {
    ss << dis(gen);
  }
  ss << "-4";
  for (i = 0; i < 3; i++) {
    ss << dis(gen);
  }
  ss << "-";
  ss << dis2(gen);
  for (i = 0; i < 3; i++) {
    ss << dis(gen);
  }
  ss << "-";
  for (i = 0; i < 12; i++) {
    ss << dis(gen);
  };
  return ss.str();
}
}


#endif //TDBSCAN_AUX_H
