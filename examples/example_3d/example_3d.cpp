//
// Created by mzoll on 01/08/2026.
//

#include <cstdlib>

#include "tdbscan_algo/common_defs.h"

#include "blib.h"
#include "connectors.h"
#include "gen_blibs.h"
#include "external/common_clib/stopwatch.h"

#include "tdbscan_algo/tdbscan_algo.h"


using namespace std;
using namespace tdbscan;
using namespace common_clib;
using namespace ex3d;

TDBScan_Algo<Blib3dWithTrace> construct_algo(const double distance_lim, const double time_lim) {
  auto limcon = new LimitingConnector(distance_lim, time_lim);

  TDBScan_Algo<Blib3dWithTrace>::TDBScan_ParameterSet params;

  params.multiplicity = 4;
  params.multiplicityTimeWindow = 2.;
  params.emergenceTimeWindow = 2.; //internally converted into SBlibWithTrace<..., tTime>::tTime::Time_t
  params.earlyMergeMultiplicityRatio = 1.; //ZeroOne: its a ratio
  params.lateMergeOverlapRatio = 1.; //ZeroOne: its a ratio

  return TDBScan_Algo<Blib3dWithTrace>(params, limcon);
}


double rand_ord() { return rand() * 100 - 50.; };
ContPos3d rand_pos() { return ContPos3d(rand_ord(), rand_ord(), rand_ord()); };
ScalarTime_t rand_time() { return ScalarTime_t(rand() % 10000); };

// std::set<Blib4d, Blib4d::TimeOrder> construct_blibs() {
std::set<Blib3dWithTrace> construct_blibs() {
  int many_blibs = 1000;

  std::set<Blib3dWithTrace> blibs;
  for (int i = 0; i < many_blibs; i++) {
    blibs.insert(Blib3dWithTrace(rand_pos(), rand_time(), Blib3dWithTrace::NOISE));
  }

  return blibs;
}


int main(int argc, char **argv) {
  auto my_algo = construct_algo({4.}, {4});

  auto blibs = construct_blibs();

  Stopwatch<std::chrono::microseconds> swatch("Process", Stopwatch<std::chrono::microseconds>::policy::start);
  auto result = my_algo.Process(blibs);
  swatch.stop();
  LOG_INFO("Processing of {} took {} {} : {} {} per blib", blibs.size(), swatch.time(), swatch.timeunitString(), (double)swatch.time()/blibs.size(), swatch.timeunitString());


  return 0;
}
