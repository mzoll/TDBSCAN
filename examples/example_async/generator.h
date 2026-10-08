//
// Created by netsu on 08/10/2026.
//

#ifndef TDBSCAN__EXASYNC__GENERATOR_H
#define TDBSCAN__EXASYNC__GENERATOR_H

#include <functional>
#include <random>
#include <utility>

#include "tdbscan_algo/common_blibs.h"

using namespace std;
using namespace tdbscan;

using namespace exasync;


namespace exasync {
class Generator {
public :
  struct generator_exhausted : std::exception {};
private:
  // create a random double on
  static double rand_double() {
    double lower_bound = 0.;
    double upper_bound = 1.;
    static std::uniform_real_distribution<double> unif(lower_bound, upper_bound);
    static std::default_random_engine re;
    return unif(re);
  }

  static std::set<ScalarBlib>
  generate_noise(const double noise_freq, const double width_fields, const double time_duration) {
    std::set<ScalarBlib> blibs;
    for (int time_step = 0; time_step < time_duration; time_step++) {
      for (int count_noise = 0; count_noise < noise_freq * width_fields; count_noise++) {
        const double pos = rand_double() * width_fields;
        const double t = rand_double() + time_step;
        blibs.insert(ScalarBlib({pos}, t));
      }
    }
    return blibs;
  }
private:
  mutable set<ScalarBlib> blib_pool_;

public:
  Generator(const double noise_freq, const double width_fields, const double time_duration)
  : blib_pool_(generate_noise(noise_freq, width_fields, time_duration)) {};
public:
  ScalarBlib next() const {
    while (!blib_pool_.empty()) {
      const auto blib = *blib_pool_.cbegin();
      blib_pool_.erase(blib_pool_.cbegin());
      return blib;
    }
    throw generator_exhausted();
  };
};

} // namespace exasync

#endif //TDBSCAN__EXASYNC__GENERATOR_H
