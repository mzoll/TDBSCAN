/**
 * \file tdbscan_algo.hh
 *
 * (c) 2026
 *
 * \author Marcel Zoll <marcel.zoll.physics@gmail.com>
 */

#pragma once

#include <algorithm>
#include <math.h>
#include <cassert>
#include <format>

#include "tdbscan_algo/dummy_logging.h"

//===========================================
//============== IMPLEMENTATION =============
//===========================================

namespace tdbscan {

//=============== namespace tdbsscan::details =================

//=============== class TDBScan_ParameterSet =================


template <class tBlib>
TDBScan_Algo<tBlib>::TDBScan_ParameterSet::TDBScan_ParameterSet() {};


//=============== class TDBScan_Algo =================================

template <class tBlib>
bool
TDBScan_Algo<tBlib>::BlibSetTimeOrder::operator()(const BlibSet &lhs, const BlibSet &rhs) const {
  return lhs.cbegin()->getTime() < rhs.cbegin()->getTime();
}

template <class tBlib>
TDBScan_Algo<tBlib>::TDBScan_Algo (
  const TDBScan_ParameterSet& params,
  Connector_t* connector) :
  params_(params),
  connector_(connector)
{
  if (params_.multiplicity<=0)
    log_fatal("Multiplicity should be greater than zero");
  if (params_.multiplicityTimeWindow<=0.0)
    log_fatal("TimeWindow should be greater than zero");
//  if (params_.acceptTimeWindow<0.0)
//    log_fatal("AcceptTimeWindow cannot be negative");
//  if (params_.rejectTimeWindow<0.0)
//    log_fatal("RejectTimeWindow cannot be negative");
//  if (params_.mergeOverlap==0)
//    log_warn("RequiredDOMOverlap configured with 0, everything will be merged");

//  if (params_.rejectTimeWindow <= params_.acceptTimeWindow)
//    log_fatal("RejectTimeWindow needs to be greater than AcceptTimeWindow");

  if (! connector_)
    log_error("No ConnectionBlock defined!");

  log_info("This is TDBScan!");
  log_debug("Leaving Init()");
};


template <class tBlib> template<class tBlibContainer>
typename TDBScan_Algo<tBlib>::BlibSetSequence
TDBScan_Algo<tBlib>::Process (const tBlibContainer& blibs) {
  log_debug("Entering Process()");
  concluded_clusters_.clear();
  active_clusters_.clear();
  emerging_clusters_.clear();

  //enforce time-order be converting to an explicit time-ordered container
  BlibSet blibs_to; // time-order
  for (const auto& b : blibs)
    blibs_to.insert(b);

  for (const auto& b : blibs_to)
    //process through machinery
    NextBlib(b);

  log_debug("Finalize");
  Finalize();

  //prepare output
  BlibSetSequence bss;
  for (const auto& c : concluded_clusters_)
    bss.insert(BlibSet(c.hits.begin(), c.hits.end() ));

  log_debug("Leaving Process()");
  return bss;
};


template <class tBlib>
void
TDBScan_Algo<tBlib>::Finalize() {
  sync_time = Time_t::max();

  emerging_clusters_.clear();

  auto ac_iter = active_clusters_.begin();
  while (ac_iter != active_clusters_.end()) {
    concluded_clusters_.push_back(*ac_iter);
    ac_iter = active_clusters_.erase(ac_iter);
  }

  //
  // //--- this is implementing the postmerge ---
  // auto ac_iter = active_clusters_.begin();
  // auto c_riter = concluded_clusters_.end();
  // while (ac_iter != active_clusters_.end()) {
  //   ac_iter->status = CausalCluster<tBlib>::CONCLUDED;
  //
  //   bool merged_any = 0;
  //   while (c_riter != concluded_clusters_.begin()) {
  //     // if we sort the list of clusters first, we could establish exit conditions faster
  //     if (c_riter->getLatestTime < ac_iter->getEarliestTime() || ac_iter->getLatestTime < c_riter->getEarliestTime()) {
  //       c_riter++;
  //       continue;
  //     }
  //     if (c_riter->nOverlap(*ac_iter)/ac_iter->count() >= params_.lateMergeOverlapRatio) {
  //       //merge
  //       c_riter->copyHits(ac_iter->getHits());
  //       merged_any = true;
  //     }
  //     ++c_riter;
  //   }
  //   if (! merged_any) {
  //     concluded_clusters_.push_back(*ac_iter);
  //   }
  //   ac_iter = active_clusters_.erase(ac_iter);
  //   ac_iter++;
  // }

  assert(active_clusters_.empty());
};


//specialize for BlibSet, which is already time-ordered
template <class tBlib>
typename TDBScan_Algo<tBlib>::BlibSetSequence
TDBScan_Algo<tBlib>::Process (const std::set<tBlib>& blibs) {
  log_debug("Entering Process()");
  emerging_clusters_.clear();
  active_clusters_.clear();
  concluded_clusters_.clear();

  //process through machinery
  for (const auto& b: blibs){
    log_debug(std::format("next Blib:: emerging: {} ; active: {} ; concluded: {}", emerging_clusters_.size(), active_clusters_.size(), concluded_clusters_.size()));
    NextBlib(b);
  }

  log_debug("Finalize");
  Finalize();

  //prepare output
  BlibSetSequence bss;
  for (const auto& c : concluded_clusters_)
    bss.insert(BlibSet(c.blibs_.begin(), c.blibs_.end() ));

  log_debug("Leaving Process()");
  return bss;
};


template <class tBlib>
bool TDBScan_Algo<tBlib>::CausallyConnected(const tBlib& b1, const tBlib& b2) const {
  return detail::CausallyConnected(b1, b2, *connector_);
};


template <class tBlib>
void TDBScan_Algo<tBlib>::NextBlib (const tBlib& b) {
  log_debug("Entering NextBlib()");
  log_trace(std::format("=== PROBE : o:{} t:{}", double(b.getOrdinate()), double(b.getTime() )));
  const auto now = b.getTime();
  // advance every each cluster in time and try to add the hit to it

  // 10. go through all active clusters  and see if blibs have fallen out of the emergence time window and multiplicity cannot be fullfilled; mark them as 'dying'; if there is nothing left mark as 'concluded'
  // 20. go through all emerging clusters and see if blibs have fallen out of the emergence time window, kill them off
  // 30. try to add the hit to existing emerging clusters; if it was added put the emerging clusters on a new_established list;
  // 40. try to add the hit to existing active clusters;
  // 50. traverse the newly established list and try to merge clusters with the active clusters
  // 9. put the hit on a newly created cluster by its own

  log_debug("Eliminating emerging clusters, adding to emerging clusters");
  // 20. go through all emerging clusters and see if blibs have fallen out of the emergence time window, if so kill the cluster off
  // 21. try to add the hit to remaining clusters; if it was added and cluster establishes (multiplicity met) put the clusters on a new_established list;
  std::list<CausalCluster<tBlib>> _newly_established_clusters;
  auto ec_iter = emerging_clusters_.begin();
  while (ec_iter != emerging_clusters_.end()) {
    const auto _n_active = ec_iter->nHitsWithinTimeWindow(  now - params_.emergenceTimeWindow, now);
    if (_n_active < ec_iter->count()) {
      log_trace("Killing one emerging cluster!");
      ec_iter = emerging_clusters_.erase(ec_iter);
      continue;
    }

    const auto success = TryInsertHit_Emergence(*ec_iter, b);
    if (success && ec_iter->count() == params_.multiplicity) {
      log_trace("Promote one emerging cluster!");
      _newly_established_clusters.push_back(*ec_iter);
      ec_iter = emerging_clusters_.erase(ec_iter);
      continue;
    }
    ec_iter++;
  }

  log_debug("Traversing active clusters");
  // 10. go through all active clusters and see if blibs have fallen out of the emergence time window and multiplicity cannot be fulfilled; mark them as 'dying'; if there is nothing left mark as 'concluded'
  auto ac_iter = active_clusters_.begin();
  while (ac_iter != active_clusters_.end()) {
    const auto _n_active = ac_iter->nHitsWithinTimeWindow(  now, Time_t::max());

    if (_n_active == 0) {
      log_trace("Active cluster concluded");
      concluded_clusters_.push_back(*ac_iter);
      ac_iter = active_clusters_.erase(ac_iter);
      continue;
    }
    else {
      log_trace("Added to active cluster");
      const auto success = TryInsertHit_Established(*ec_iter, b);
    }
    ac_iter++;
  }


  log_debug("Early merge");
  // 50. traverse the newly established list and try to merge clusters with the active clusters
  auto nec_iter = _newly_established_clusters.begin();
  ac_iter = active_clusters_.begin();
  while (nec_iter != _newly_established_clusters.end()) {
    while (ac_iter != active_clusters_.end()) {

      //TODO implement the early merge criteria

      if (nec_iter->isSubsetOf(*ac_iter)) {
        nec_iter = _newly_established_clusters.erase(nec_iter);
        ac_iter = active_clusters_.begin();
        continue;
      }
      ac_iter++;
    }
    //established cluster is a genuinely new cluster
    nec_iter++;
    ac_iter = active_clusters_.begin();
  }
  active_clusters_.insert(active_clusters_.end(), _newly_established_clusters.begin(), _newly_established_clusters.end());
  _newly_established_clusters.clear();

  log_debug("Create self-contained cluster");
  // 9. put the hit on a newly created cluster by its own
  emerging_clusters_.insert(emerging_clusters_.end(), CausalCluster(b));

  log_debug("Leaving NextHit()");
}


template <class tBlib>
bool TDBScan_Algo<tBlib>::TryInsertHit_Emergence(
  CausalCluster<tBlib>& c,
  const tBlib& b) {
  log_debug("Entering TryInsertHit_Emergence()");
  //all blibs in the emerging Cluster must connect

  for (const auto& cb : c.blibs_) {
    log_trace("loop");

    if (cb.timeDiff(b) >= params_.emergenceTimeWindow) {
      /// we are past the timeframe;
      return false;
    }
    if (not CausallyConnected(cb, b))
      return false;
  }
  log_trace("Adding");
  c.blibs_.insert(c.blibs_.end(), b);

  log_debug("Leaving TryInsertHit_Established()");
  return true;
}


template <class tBlib>
bool TDBScan_Algo<tBlib>::TryInsertHit_Established(
  CausalCluster<tBlib>& c,
  const tBlib& b) {
  log_debug("Entering TryInsertHit_Established()");
  auto _result = false;

  int active_connectees = 0;
  auto cb_iter = c.blibs_.cend();
  while (cb_iter != c.blibs_.cbegin()) {
    if (cb_iter->timeDiff(b) >= params_.multiplicityTimeWindow) {
      /// we are past the timeframe;
      break;
    }
    active_connectees += CausallyConnected(*cb_iter, b);
  }
  if (active_connectees >= params_.multiplicity) {
    c.blibs_.insert(c.blibs_.end(), b);
    _result = true;
  }

  log_debug("Leaving TryInsertHit_Established()");
  return _result;
};



} //namespace tdbscan
