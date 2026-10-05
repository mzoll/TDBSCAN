//
// Created by netsu on 05/10/2026.
//

#pragma once

#include "tdbscan_async.h"


template<class tBlib>
TDBScan_Algo<tBlib>::BlibSet
TDBScan_Algo<tBlib>::ObtainOutput() {
  sig_output_available_.wait(output_access_lock_);
  const auto output= output_queue_.pop();
  output_access_lock_.unlock();
  return output.getBlibs();
}



template<class tBlib>
void
TDBScan_Algo<tBlib>::FeedBlib(const tBlib &b) {
  if (finalize_requested_)
    throw std::logic_error("::FeedBlib() requested after ::Finalize() has been called;");

  mtx_input_access_.lock();
  input_queue_.push(b);
  mtx_input_access_.unlock();
  sig_input_available_.notify_all();
}

template<class tBlib>
void
TDBScan_Algo<tBlib>::Finalize() {
  finalize_requested_ = true;


  //needs to close the input queue
  //wait: let the rest of input blibs be processed
  //
  sig_ready_for_finalize_.wait(glob_seqaccess_lock_);

  if (!concluded_clusters_.empty) {
    //transfer all concluded clusters to the output queue
    output_access_lock_.lock();
    while (!concluded_clusters_.empty()) {
      output_queue_.push(concluded_clusters_.pop_front());
    }
    output_access_lock_.unlock();
    sig_output_available_.notify_all();
  }

  assert(active_clusters_.empty());
};
