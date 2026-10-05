//
// Created by netsu on 05/10/2026.
//

#ifndef TDBSCAN_TDBSCAN_ASYNC_H
#define TDBSCAN_TDBSCAN_ASYNC_H

#include <thread>

#include "tdbscan_algo.h"

#include "external/common_clib/InterruptableQueue.hpp"

namespace tdbscan {
/**
 * The main algorithm class
 * Needs to be configured with a connector and a parameter set then can be fed with a sequence of hits
 */
template<class tBlib>
class TDBScan_AsyncMachine : protected TDBScan_Algo<tBlib> {


  ///an FIFO for putting blibs on
  common_clib::threading::InterruptableQueue<tBlib> input_queue_;
  ///an FIFO for laying the output on
  common_clib::threading::InterruptableQueue<typename TDBScan_Algo<tBlib>::CausalCluster_t> output_queue_;

  std::mutex mtx_output_access_;
  /// coordinates access to the `output_queue_`
  std::unique_lock<std::mutex> output_access_lock_{mtx_output_access_};

  /// signals the availability of at least one available output
  std::condition_variable sig_input_available_;
  std::condition_variable sig_output_available_;
  std::condition_variable sig_ready_for_finalize_;

  std::thread* driving_tread_{nullptr};

  bool finalize_requested_{false};

  struct Telemetry {
    unsigned int n_blibs_{0};
    unsigned int n_output_clusters_{0};

  } telemetry_;

public: // --- THE REAL MACHINERY ---
  //===================
  // Internal Methods
  //===================

  /** Asynchronously put one blib on the input queue
   * @param b the blib to feed
   */
  void FeedBlib(const tBlib &b);

  void Crank();

  /**
 * Obtain one output cluster as it becomes available
 * @return a set of blibs which is a cluster
 */
  BlibSet ObtainOutput();


  void Finalize();

public:
  //machine control
  ///start the machine, aka the internal thread
  void start() noexcept;
  ///stop the internal thread
  void stop() noexcept;

public:
  /**
 * Constructor from a ParameterSet and Connector
 * @param params Collection of Parameters which govern the algorithm
 * @param connector Pointer to a Connector, which facilítates the comparison of hits
 */
  TDBScan_AsyncMachine(
    const TDBScan_ParameterSet &params,
    const Connector_t *connector);

  /**
   * Destructor
   * hold the internal thread, when started
   *
   */
  ~TDBScan_AsyncMachine();
};
} //namespace tdbscan;

#include "tdbscan_async.hh"

#endif //TDBSCAN_TDBSCAN_ASYNC_H
