//
// Created by netsu on 05/10/2026.
//

#ifndef TDBSCAN_TDBSCAN_ASYNC_H
#define TDBSCAN_TDBSCAN_ASYNC_H

#include "tdbscan_algo.h"


/**
 * The main algorithm class
 * Needs to be configured with a connector and a parameter set then can be fed with a sequence of hits
 */
template<class tBlib>
class TDBScan_AsyncMachine : protected TDBScan_Algo {


  ///an FIFO for putting blibs on
  std::queue<CausalCluster_t> input_queue_;
  ///an FIFO for laying the output on
  std::queue<CausalCluster_t> output_queue_;


  tdbscan::TDBScan_Algo<> algo_;



  std::mutex mtx_output_access_;
  /// coordinates access to the `output_queue_`
  std::unique_lock<std::mutex> output_access_lock_{mtx_output_access_};

  /// signals the availability of at least one available output
  std::condition_variable sig_input_available_;
  std::condition_variable sig_output_available_;
  std::condition_variable sig_ready_for_finalize_;



public: // --- THE REAL MACHINERY ---
  //===================
  // Internal Methods
  //===================

  /** Asynchronously put one blib on the input queue
   * @param b the blib to feed
   */
  void FeedBlib(const tBlib &b);

  /**
 * Obtain one output cluster as it becomes available
 * @return a set of blibs which is a cluster
 */
  BlibSet ObtainOutput();


  void Finalize();


public:
  /**
 * Constructor from a ParameterSet and Connector
 * @param params Collection of Parameters which govern the algorithm
 * @param connector Pointer to a Connector, which facilítates the comparison of hits
 */
  TDBScan_AsyncMachine(
    const TDBScan_ParameterSet &params,
    const Connector_t *connector) : TDBScan_Algo(params, connector) {};
};


#include "tdbscan_async.hh"

#endif //TDBSCAN_TDBSCAN_ASYNC_H
