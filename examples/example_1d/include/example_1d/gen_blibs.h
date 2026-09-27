//
// Created by netsu on 27/09/2026.
//

#ifndef TDBSCAN__EX1D__GEN_BLIBS_H
#define TDBSCAN__EX1D__GEN_BLIBS_H

#include <random>
#include "blib.h"

#include "tdbscan_algo/auxilary/trivial_logging.h"

// create a random double on
double rand_double() {
  double lower_bound = 0.;
  double upper_bound = 1.;
  static std::uniform_real_distribution<double> unif(lower_bound,upper_bound);
  static std::default_random_engine re;
  return unif(re);
}

namespace ex1d {
using namespace tdbscan;

/* ========================= Create the scenario =======================
 * For a demonstration and testing scenario Blibs need to be generated which stem from two distinguishable sources:
 * Signal and Noise.
 * In Our case of a 1d- grid, the signal will be a represented by a box of finite size moving sideways
 * with a given inertia generating blibs according to a brightness. Noise on the other hand will be random noise
 * uniformly generated on the fields of that grid.
 */
std::set<SBlibWithTrace>
generate_noise(const double noise_freq, const double width_fields, const double time_duration) {
  std::set<SBlibWithTrace> blibs;
  for (int time_step = 0; time_step < time_duration; time_step++) {
    for (int count_noise = 0; count_noise < noise_freq * width_fields; count_noise++) {
      const double pos = rand_double() * width_fields;
      const double t = rand_double() + time_step;
      blibs.insert(SBlibWithTrace({pos}, t, SBlibWithTrace::NOISE));
    }
  }
  return blibs;
}


/**
 * generate blibs from a box moving over a 1d space left to right
 *
 * example: box_size: 4, inertia: 2, ledge_start_pos: 0, brightness: 0.5, time_duration 5
 * --[xx  ]-------------------
 * --.--[ x x]----------------
 * --.----[x  x]--------------
 * --.------[ xx ]------------
 * --.--------[  xx]----------
 *
 * The box will be reflected on the right and left side
 *
 * @param box_size
 * @param inerta
 * @param start_pos,
 * @param brightness a measure of the signal frequency in one unit-volume of the box
 * @param field_size,
 * @param time_duration
 * @return
 */
std::set<SBlibWithTrace>
generate_moving_box(
  const double box_size,
  const double inertia,
  const double start_pos, // =0.,
  const double brightness, // =1.,
  const double field_size,
  const double time_duration) {
  std::set<SBlibWithTrace> blibs;

  const int n_blibs = static_cast<int>(time_duration * box_size * brightness);

  for (int i_blib = 0; i_blib < n_blibs; i_blib++) {
    const auto t = rand_double() * time_duration;
    const auto center_pos_ind = (t*inertia + start_pos) / field_size;

    const double center_pos = field_size * (int(center_pos_ind) % 2 == 0 ?  center_pos_ind - std::floor(center_pos_ind) : 1. - (center_pos_ind - std::floor(center_pos_ind)));

    const double blib_pos_inbox = box_size * (rand_double() - 1/2.);

    const auto blib_pos = center_pos + blib_pos_inbox;

    if (blib_pos<0. || blib_pos > field_size)
      continue;

    blibs.insert(SBlibWithTrace({blib_pos}, t, SBlibWithTrace::SIGNAL));
  }

  return blibs;
}


/**
 * Generate Blibs for our scenario
 *
 * A bright box moves left to right in
 * @param time_duration
 * @param width_fields
 * @param brightness
 * @param noise_contamination
 * @return
 */
std::set<SBlibWithTrace>
gernerate_blibs( const double time_duration=50, const double width_fields = 100, const double brightness= 1., const double noise_contamination = 0.1) {
  std::set<SBlibWithTrace> blibs;

  const auto _box_blibs = generate_moving_box(5., 1., 0., brightness, width_fields, time_duration);
  LOG_INFO("Generated {} BOX blibs", _box_blibs.size());
  blibs.insert(_box_blibs.cbegin(), _box_blibs.cend());

  const auto _noise_blibs = generate_noise(noise_contamination*brightness, width_fields, time_duration);
  LOG_INFO("Generated {} NOISE blibs", _noise_blibs.size());
  //blibs.insert(_noise_blibs.cbegin(), _noise_blibs.cend());
  return blibs;
}
} //namespace ex1d

#endif //TDBSCAN__EX1D__GEN_BLIBS_H
