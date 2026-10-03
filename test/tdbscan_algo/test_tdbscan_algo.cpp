//
// Created by netsu on 03/10/2026.
//

#include "gtest/gtest.h"
#include "tdbscan_algo/common_blibs.h"
#include "tdbscan_algo/tdbscan_algo.h"
#include <random>

using namespace tdbscan;


// create a random double on
double rand_double() {
  double lower_bound = 0.;
  double upper_bound = 1.;
  static std::uniform_real_distribution<double> unif(lower_bound, upper_bound);
  static std::default_random_engine re;
  return unif(re);
}

/* ========================= Create the scenario =======================
 * For a demonstration and testing scenario Blibs need to be generated which stem from two distinguishable sources:
 * Signal and Noise.
 * In Our case of a 1d- grid, the signal will be a represented by a box of finite size moving sideways
 * with a given inertia generating blibs according to a brightness. Noise on the other hand will be random noise
 * uniformly generated on the fields of that grid.
 */
std::set<ScalarBlib>
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





class Connector_TRUE : public ConnectorSingle<ScalarBlib> {
public:
  bool eval(const ScalarBlib &h1, const ScalarBlib &h2) const override { return true; };
  Connector_TRUE() : ConnectorSingle("True") {};
};

class Connector_FALSE : public ConnectorSingle<ScalarBlib> {
public:
  bool eval(const ScalarBlib &h1, const ScalarBlib &h2) const override { return false; };
  Connector_FALSE() : ConnectorSingle("False") {};
};


/* ============================ Constructing the TDBscan algorithm instance ==================
 * the main algorithm is constructed by creating instances of the limiters defined in the previous set,
 * and providing a set of Algorithm parameters
 */
TDBScan_Algo<ScalarBlib> construct_algo(const bool everything_connected) {

  auto con = everything_connected ?
    dynamic_cast<ConnectorSingle<ScalarBlib>*>( new Connector_TRUE) :
    dynamic_cast<ConnectorSingle<ScalarBlib>*>(new Connector_FALSE);

  TDBScan_Algo<ScalarBlib>::TDBScan_ParameterSet params;

  params.multiplicity = 4;
  params.multiplicityTimeWindow = 2.; //internally converted into SBlibWithTrace<..., tTime>::tTime::Time_t
  params.emergenceTimeWindow = 2.; //internally converted into SBlibWithTrace<..., tTime>::tTime::Time_t
  params.earlyMergeMultiplicityRatio = 1.; //ZeroOne: its a ratio
  params.lateMergeOverlapRatio = 1.; //ZeroOne: its a ratio

  return TDBScan_Algo(params, con);
}



TEST(AlgoTest, HappyTest_TRUE) {
  const auto blibs = generate_noise(1, 10, 100);

  auto my_algo = construct_algo(true);

  const auto result = my_algo.Process(blibs);
}

TEST(AlgoTest, HappyTest_FALSE) {
  const auto blibs = generate_noise(1, 10, 100);

  auto my_algo = construct_algo(true);

  const auto result = my_algo.Process(blibs);
}
