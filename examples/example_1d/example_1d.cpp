//
// Created by mzoll on 01/08/2026.
//

#include "blib.h"
#include "connectors.h"
#include "gen_blibs.h"
#include "helpers.h"

#include "tdbscan_algo/tdbscan_algo.h"

using namespace std;
using namespace tdbscan;
using namespace ex1d;


/* ============================ Constructing the TDBscan algorithm instance ==================
 * the main algorithm is constructed by creating instances of the limiters defined in the previous set,
 * and providing a set of Algorithm parameters
 */
TDBScan_Algo<SBlibWithTrace> construct_algo(const double distance_lim, const double time_lim) {
  auto limcon = new LimitingConnector(distance_lim, time_lim);

  TDBScan_Algo<SBlibWithTrace>::TDBScan_ParameterSet params;

  params.multiplicity = 4;
  params.multiplicityTimeWindow = 2.; //internally converted into SBlibWithTrace<..., tTime>::tTime::Time_t
  params.emergenceTimeWindow = 2.; //internally converted into SBlibWithTrace<..., tTime>::tTime::Time_t
  params.earlyMergeMultiplicityRatio = 1.; //ZeroOne: its a ratio
  params.lateMergeOverlapRatio = 1.; //ZeroOne: its a ratio

  return TDBScan_Algo(params, limcon);
}


int main(int argc, char **argv) {
  auto my_algo = construct_algo(2., 0.5);

  LOG_INFO("Generate blibs");
  const auto blibs = gernerate_blibs(50, 100, 3, 0.1);

  LOG_INFO("Processing nBlibs: {} (purity {:.3f})", blibs.size(), calculate_signal_purity(blibs));
  //take first 3
  std::set<SBlibWithTrace> _blibs;
  auto iter = blibs.begin();
  for (int i = 0; i < 100; i++) {
    _blibs.insert(*iter);
    LOG_TRACE("Sample : {}", *iter);
    ++iter;
  }


  const auto result = my_algo.Process(_blibs);

  LOG_INFO("Generated nClusters: {}", result.size());

  auto r_citer = result.cbegin();
  for (int i = 0; i < 5; i++) {
    if (r_citer == result.cend())
      break;
    LOG_INFO("Cluster {} size: {} (purity {:.3f})", i, r_citer->size(), calculate_signal_purity(*r_citer));
    for (const auto &b: *r_citer) {
      LOG_INFO("das {}", b)
    }

    r_citer++;
  }

  // for (const auto& c : result) {
  // 	LOG_INFO(std::format("Size: {}", c.size()));
  // }
}
