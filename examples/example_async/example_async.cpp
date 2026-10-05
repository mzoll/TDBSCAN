//
// Created by mzoll on 01/08/2026.
//

#include "blib.h"
#include "connectors.h"
#include "gen_blibs.h"
#include "helpers.h"
#include "external/common_clib/interrupt.h"

#include "tdbscan_algo/tdbscan_algo.h"
#include "external/common_clib/stopwatch.h"
#include "tdbscan_algo/common_blibs.h"
#include "tdbscan_algo/tdbscan_async.h"

using namespace std;
using namespace tdbscan;
using namespace common_clib;
using namespace ex1d;

using namespace std::chrono_literals;

/* ============================ Constructing the TDBscan algorithm instance ==================
 * the main algorithm is constructed by creating instances of the limiters defined in the previous set,
 * and providing a set of Algorithm parameters
 */
TDBScan_AsyncMachine<ScalarBlib> construct_algo(const double distance_lim, const double time_lim) {
  auto limcon = new LimitingConnector(distance_lim, time_lim);

  TDBScan_Algo<ScalarBlib>::TDBScan_ParameterSet params;

  params.multiplicity = 4;
  params.multiplicityTimeWindow = 2.; //internally converted into SBlibWithTrace<..., tTime>::tTime::Time_t
  params.emergenceTimeWindow = 2.; //internally converted into SBlibWithTrace<..., tTime>::tTime::Time_t
  params.earlyMergeMultiplicityRatio = 1.; //ZeroOne: its a ratio
  params.lateMergeOverlapRatio = 1.; //ZeroOne: its a ratio

  return TDBScan_AsyncMachine<ScalarBlib>(params, limcon);
}


struct generator_exhausted : std::exception {};

int main(int argc, char **argv) {
  auto async_algo = construct_algo(2., 0.5);

  LOG_INFO("Generate blibs");
  std::list<ScalarBlib> blibs = generate_blibs(5000, 100, 3, 0.1);





  auto obfuscated_blib_generator = [&blibs]() {
    if (blibs.empty())
      throw generator_exhausted();
    const auto b = blibs.front();
    blibs.pop_front();
    std::this_thread::sleep_for(50ms);
    return b;
  };

  std::list<std::set<ScalarBlib>> clusters;

  bool end_feed = false;

  auto feed = [obfuscated_blib_generator, &async_algo, &end_feed]() {
    while (!blibs.empty() and !end_feed) {
      try {
        async_algo.FeedBlib(obfuscated_blib_generator());
      } catch (generator_exhausted) {
        break;
      }
    }
    end_feed = true;
  };

  auto consume = [&clusters, &async_algo]() {
    while (true)
      try {
        clusters.push_back(async_algo.ObtainCluster());
      } catch (common_clib::threadsafe::interrupt_exception) {
        // output consumption was interrupted
        break;
      }
  };

  auto alert = nullptr; //implement this

  // start the machine
  async_algo.start();

  std::thread feeder(feed);
  std::thread consumer(consume);

  std::this_thread::sleep_for(1s);

  async_algo.stop(); // waiting consumers on async_algo::ObtainCluster() are notified and are resurfacing.
  consumer.join(); // consumer has ended execution, lasso in the thread

  end_feed = true; // causes the feeder to end execution on next loop iteration, if not already exhausted

  feeder.join(); // feeder has ended execution, lasso in the thread

  async_algo.Finalize();  // sequentially work off, what is

  while (async_algo.more_output())
    clusters.push_back(async_algo.ObtainCluster());

}
