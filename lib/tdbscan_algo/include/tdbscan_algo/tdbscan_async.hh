//
// Created by netsu on 05/10/2026.
//

#pragma once

#include "tdbscan_async.h"

namespace tdbscan {

template<class tBlib>
TDBScan_AsyncMachine<tBlib>::TDBScan_AsyncMachine(
  const TDBScan_Algo<tBlib>::TDBScan_ParameterSet &params,
  const TDBScan_Algo<tBlib>::Connector_t *connector) : algo_(params, connector) {};

template<class tBlib>
TDBScan_AsyncMachine<tBlib>::~TDBScan_AsyncMachine() {
  stop();
};

template<class tBlib>
void
TDBScan_AsyncMachine<tBlib>::start() noexcept {
  if (input_queue_.is_closed())
    input_queue_.reopen();
  if (!driving_tread_)
    driving_tread_ = new thread(&TDBScan_AsyncMachine::Crank, this);
}

template<class tBlib>
void
TDBScan_AsyncMachine<tBlib>::stop() noexcept {
  LOG_DEBUG("Stop called for TDBSCAN Machine");
  if (driving_tread_) {
    input_queue_.close();
    driving_tread_->join();
    input_queue_.reopen();
    delete driving_tread_;
    driving_tread_ = nullptr;
  }
}

template<class tBlib>
void
TDBScan_AsyncMachine<tBlib>::FeedBlib(const tBlib &b) {
  if (finalize_requested_)
    throw std::logic_error("::FeedBlib() requested after ::Finalize() has been called;");
  input_queue_.push(b);
}

template<class tBlib>
void TDBScan_AsyncMachine<tBlib>::Crank() {
  while (true) {
    try {
      //within here there is serialized access
      algo_.NextBlib(input_queue_.pop());  //this blocks until a blib actually becomes available on the queue
      ++telemetry_.n_blibs_;

      while (!algo_.concluded_clusters_.empty()) {
        output_queue_.push(algo_.concluded_clusters_.front());
        algo_.concluded_clusters_.pop_front();
        ++telemetry_.n_output_clusters_;
      }

    } catch (typename common_clib::threading::InterruptableQueue<tBlib>::InterruptIsSet& i) {
      // consumption has been interrupted
      LOG_ERROR("Input queue is blocked");
      break;
    } catch (...) {
      // consumption has been interrupted by something unforseen
      throw;
    }
  }
  LOG_ERROR("Ending driving thread execution");
}


template<class tBlib>
TDBScan_Algo<tBlib>::BlibSet
TDBScan_AsyncMachine<tBlib>::ObtainCluster() {
  return output_queue_.pop().getBlibs();
}

template<class tBlib>
void
TDBScan_AsyncMachine<tBlib>::Finalize(const bool consume_inputs) {
  finalize_requested_ = true;
  //bleed the input-queue dry

  input_queue_.block(); //makes the driving thread resurface and stop

  if (consume_inputs) {
    for (auto& b : input_queue_.exhaust()) {
      NextBlib(input_queue_.pop());  //this blocks until a blib actually becomes available on the queue
      ++telemetry_.n_blibs_;
    }
  }
  TDBScan_Algo<tBlib>::Finalize();

  while (!algo_.concluded_clusters_.empty()) {
    output_queue_.push(algo_.concluded_clusters_.pop_front());
    ++telemetry_.n_output_clusters_;
  }
}

} //namespace tdbscan;
