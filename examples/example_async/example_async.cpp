//
// Created by mzoll on 01/08/2026.
//


#include "connector.h"
#include "generator.h"
#include "agents.h"

#include "tdbscan_algo/tdbscan_async.h"

using namespace std;
using namespace tdbscan;

using namespace std::chrono_literals;

using namespace exasync;



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



int main(int argc, char **argv) {
  auto machine = construct_algo(2., 0.5);

  LOG_INFO("Generate blibs");

  Generator gen(20., 10, 500);

  struct feeder_task {
    const Generator* gen_ptr_;
    TDBScan_AsyncMachine<ScalarBlib>* const machine_ptr_;
    feeder_task(const Generator* gen_ptr, TDBScan_AsyncMachine<ScalarBlib>* machine_ptr)
    : gen_ptr_(gen_ptr), machine_ptr_(machine_ptr) {};
    void operator()() const {machine_ptr_->FeedBlib(gen_ptr_->next());}
  } _feeder_task(&gen, &machine);

  Feeder feeder(_feeder_task);

  auto consumer_task = [&machine]() {
    machine.ObtainCluster();
  };

  Consumer consumerA(consumer_task), consumerB(consumer_task);

  //starting all parts of the machine
  machine.start();
  //consumerA.start();
  //consumerB.start();
  feeder.start();

  std::this_thread::sleep_for(1s);
  machine.stop();
  std::this_thread::sleep_for(1s);
  feeder.stop();
  //consumerA.stop(); //does not surface!
  //consumerB.stop(); //does not surface!

  //shutdown of all parts by destructor
}
