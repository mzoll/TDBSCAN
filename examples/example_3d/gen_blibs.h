//
// Created by netsu on 27/09/2026.
//

#ifndef TDBSCAN__EXAMPLE_3D__GEN_BLIBS_H
#define TDBSCAN__EXAMPLE_3D__GEN_BLIBS_H

#include "blib.h"
#include <random>

// create a random double on
double rand_double() {
  static std::uniform_real_distribution<double> unif(0., 1.);
  static std::default_random_engine re;
  return unif(re);
}

namespace ex3d {
/**
 * generate blibs from a box that makes random appears at random positions
 * @note number of realised appearances might not be the number requested. But is guarantied in the limit of great numbers.
 *
 * example: box_size: 4, brightness: 0.5, time_duration 5, field_size 20
 * --xxoo-------------------
 * --oxox-------------------
 * ---------xoox------------
 * ---------oxxo------------
 * ----------------ooxx-----
 *
 * The box will be reflected on the right and left side
 *
 * @param box_size
 * @param brightness a measure of the signal frequency in one unit-volume of the box
 * @param appearances number of appearances to make
 * @param field_size,
 * @param time_duration
 * @return blibs generated
 */
std::set<Blib3dWithTrace>
generate_appearing_box(
  const double box_size,
  const double brightness, // =1.,
  const int appearances,
  const double field_size,
  const double time_duration) {
  std::set<Blib3dWithTrace> blibs;

  static std::default_random_engine re;
  const double mu = appearances/time_duration;
  std::normal_distribution<double> gauss(mu, sqrt(mu));

  double t_idx = 0.;
  while (t_idx < time_duration) {
    const auto a_duration = gauss(re);
    const auto a_center_pos = rand_double() * field_size;

    const int n_blibs = static_cast<int>(a_duration * box_size * brightness);
    for (int i=0; i<n_blibs; i++) {
      const auto t = rand_double() * a_duration + t_idx;
      if (t > time_duration)
        continue;

      auto rand_pos = [box_size, a_center_pos](){return box_size * (rand_double() - 1 / 2.) + a_center_pos;};
      auto x = box_size * (rand_double() - 1 / 2.) + a_center_pos;;
      auto y = box_size * (rand_double() - 1 / 2.) + a_center_pos;;
      auto z = box_size * (rand_double() - 1 / 2.) + a_center_pos;;

      if ((0. > x || x < field_size) || (0. > y || y < field_size) || (0. > z || z < field_size))
        continue;

      blibs.insert(Blib3dWithTrace({x,y,z}, {t}, Blib3dWithTrace::SIGNAL));
    }
    t_idx += a_duration;
  }
  return blibs;
};
} //namespace ex3d

#endif //TDBSCAN__EXAMPLE_3D__GEN_BLIBS_H
