/**
 * \file tdbscan_algo.hh
 *
 * (c) 2026
 *
 * \author Marcel Zoll <marcel.zoll.physics@gmail.com>
 */

#pragma once

#include <algorithm>
#include <cmath>
#include <cassert>
#include <format>


#include "auxilary/trivial_logging.h"

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
  const Connector_t* const connector) :
  params_(params),
  connector_(connector)
{
  if (params_.multiplicity<=0)
    LOG_FATAL("Multiplicity should be greater than zero");
  if (params_.multiplicityTimeWindow<=0.0)
    LOG_FATAL("TimeWindow should be greater than zero");
//  if (params_.acceptTimeWindow<0.0)
//    LOG_FATAL("AcceptTimeWindow cannot be negative");
//  if (params_.rejectTimeWindow<0.0)
//    LOG_FATAL("RejectTimeWindow cannot be negative");
//  if (params_.mergeOverlap==0)
//    LOG_WARNn("RequiredDOMOverlap configured with 0, everything will be merged");

//  if (params_.rejectTimeWindow <= params_.acceptTimeWindow)
//    LOG_FATAL("RejectTimeWindow needs to be greater than AcceptTimeWindow");

  if (! connector_)
    LOG_ERROR("No ConnectionBlock defined!");

  LOG_INFO("This is TDBScan!");
  LOG_DEBUG("Leaving Init()");
};


template <class tBlib> template<class tBlibContainer>
typename TDBScan_Algo<tBlib>::BlibSetSequence
TDBScan_Algo<tBlib>::Process (const tBlibContainer& blibs) {
  LOG_DEBUG("Entering Process()");
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

  LOG_DEBUG("Finalize");
  Finalize();

  //prepare output
  BlibSetSequence bss;
  for (const auto& c : concluded_clusters_)
    bss.insert(BlibSet(c.hits.begin(), c.hits.end() ));

  LOG_DEBUG("Leaving Process()");
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
  LOG_DEBUG("Entering Process()");
  emerging_clusters_.clear();
  active_clusters_.clear();
  concluded_clusters_.clear();

  //process through machinery
  for (const auto& b: blibs){
    NextBlib(b);
  }

  LOG_DEBUG("Finalize");
  Finalize();

  //prepare output
  BlibSetSequence bss;
  for (const auto& c : concluded_clusters_)
    bss.insert(BlibSet(c.blibs_.begin(), c.blibs_.end() ));

  LOG_DEBUG("Leaving Process()");
  return bss;
};


template <class tBlib>
bool TDBScan_Algo<tBlib>::CausallyConnected(const tBlib& b1, const tBlib& b2) const {
  return detail::CausallyConnected(*connector_, b1, b2);
};


template <class tBlib>
void TDBScan_Algo<tBlib>::NextBlib (const tBlib& b) {
  LOG_DEBUG("Entering NextBlib()");
  LOG_TRACE(">>> NEXT BLIB : {}", b );
  LOG_TRACE("current:: emerging: {} ; active: {} ; concluded: {}", emerging_clusters_.size(), active_clusters_.size(), concluded_clusters_.size());
  const auto now = b.getTime();
  // advance every each cluster in time and try to add the hit to it

  // 10. go through all emerging clusters and see if blibs have fallen out of the emergence time window, kill them off;
  // 20. go through all active clusters and see if blibs have fallen out of the time window and multiplicity cannot be fulfilled; if there is nothing left mark as 'concluded'
  // 30. try to add the hit to existing emerging clusters; if it was added, put the emerging clusters that fulfill multiplicity on a 'new_established' list;
  // 40. try to add the hit to existing active clusters;
  // 50. traverse the newly established list and try to merge clusters with the active clusters
  // 9. put the hit on a newly created cluster by its own

  LOG_DEBUG("Eliminating emerging clusters, adding to emerging clusters");
  // 10. go through all emerging clusters and see if blibs have fallen out of the emergence time window, if so kill the cluster off
  std::list<CausalCluster<tBlib>> _newly_established_clusters;
  auto ec_iter = emerging_clusters_.begin();
  while (ec_iter != emerging_clusters_.end()) {
    const auto _n_active = ec_iter->nHitsWithinTimeWindow(  now - params_.emergenceTimeWindow, now);
    if (_n_active < ec_iter->count()) {
      LOG_TRACE("Killing off one emerging cluster!");
      ec_iter = emerging_clusters_.erase(ec_iter);
      continue;
    }

    // 11. try to add the hit to remaining clusters; if it was added and cluster establishes (multiplicity met) put the clusters on a new_established list;
    const auto success = TryInsertHit_Emergence(*ec_iter, b);
    if (success && ec_iter->count() == params_.multiplicity) {
      LOG_TRACE("Promote one emerging cluster!");
      _newly_established_clusters.push_back(*ec_iter);
      ec_iter = emerging_clusters_.erase(ec_iter);
      continue;
    }
    ++ec_iter;
  }

  LOG_DEBUG("Create self-contained cluster");
  // 12. put the hit on a newly created emerging cluster by its own
  emerging_clusters_.insert(emerging_clusters_.end(), CausalCluster(b));


  LOG_DEBUG("Traversing active clusters");
  // 20. go through all active clusters and see if blibs have fallen out of the time window and multiplicity cannot be fulfilled; if there is nothing left, move to 'concluded'
  // 21. try to add hit to any Active cluster that has sufficient causal evidence
  auto ac_iter = active_clusters_.begin();
  while (ac_iter != active_clusters_.end()) {
    const auto _n_active = ac_iter->nHitsWithinTimeWindow(  now-params_.emergenceTimeWindow, now);
    if (_n_active == 0) {
      LOG_TRACE("Active cluster concluded");

      // TODO late merge mechanism here
      LOG_DEBUG("Late merge mechanism");
      auto acother_iter = active_clusters_.begin();
      bool was_reabsorbed = false;
      while (acother_iter != active_clusters_.end()) {
        if (ac_iter == acother_iter) {
          // its the cluster itself
          ++acother_iter;
          continue;
        }

        if (ac_iter->nOverlap(*acother_iter)/ std::min(ac_iter->count(), acother_iter->count()) >= params_.lateMergeOverlapRatio) {
          LOG_TRACE("this cluster can be reabsorbed ...");
          was_reabsorbed = true;
          acother_iter->copyBlibs(*ac_iter);
        }
        ++acother_iter;
      }
      if (was_reabsorbed) {
        LOG_TRACE("... and was thus erased")
        ac_iter = active_clusters_.erase(ac_iter);
        continue;
      }

      LOG_TRACE("... and was pushed to the concluded clusters");
      concluded_clusters_.push_back(*ac_iter);
      ac_iter = active_clusters_.erase(ac_iter);
      continue;
    }

    LOG_TRACE("Try Adding Blib to active cluster");
    const auto success = TryInsertHit_Established(*ac_iter, b);
    ++ac_iter;
  }

  LOG_DEBUG("Early merge: iterating _newly_established_clusters: {}", _newly_established_clusters.size());
  LOG_TRACE("current:: emerging: {} ; active: {} ; concluded: {}", emerging_clusters_.size(), active_clusters_.size(), concluded_clusters_.size());
  // 50. traverse the newly established list and try to merge clusters with the active clusters
  auto nec_citer = _newly_established_clusters.cbegin();
  ac_iter = active_clusters_.begin();
  while (nec_citer != _newly_established_clusters.cend()) {
    while (ac_iter != active_clusters_.end() && nec_citer != _newly_established_clusters.cend()) {

      //TODO implement the early merge criteria

      // if there is overlap except in one blib, the one holdout blib might be noise hit that was just , but the
      const auto overlap = nec_citer->nOverlap(*ac_iter);
      if ( overlap >= params_.multiplicity - params_.earlyMergeRejectionHoldout ) {
        LOG_TRACE("this is a subset;");
        nec_citer = _newly_established_clusters.erase(nec_citer);
        ac_iter = active_clusters_.begin();
        continue;
      }
      LOG_TRACE("this is NOT a subset;");
      ++ac_iter;
    }
    //established cluster is a genuinely new cluster
    ++nec_citer;
    ac_iter = active_clusters_.begin();
  }
  if (!_newly_established_clusters.empty()) {
    LOG_DEBUG("Inserting {} genuine new active clusters", _newly_established_clusters.size());
    active_clusters_.insert(active_clusters_.end(), _newly_established_clusters.begin(), _newly_established_clusters.end());
    _newly_established_clusters.clear();
  }

  LOG_TRACE("current:: emerging: {} ; active: {} ; concluded: {}", emerging_clusters_.size(), active_clusters_.size(), concluded_clusters_.size());
  LOG_DEBUG("Leaving NextHit()");
}


template <class tBlib>
bool TDBScan_Algo<tBlib>::TryInsertHit_Emergence(
  CausalCluster<tBlib>& c,
  const tBlib& b) {
  LOG_DEBUG("Entering TryInsertHit_Emergence()");
  //all blibs in the emerging Cluster must connect

  for (const auto& cb : c.blibs_) {
    const auto timediff = cb.timeTo(b);
    if (timediff>= params_.emergenceTimeWindow ) {
      LOG_TRACE("Timediff past the allowed emergenceTimeWindow({}): {}", double(params_.emergenceTimeWindow), double(timediff));
      /// we are past the timeframe;
      return false;
    }
    if (not CausallyConnected(cb, b)) {
      LOG_TRACE("Blibs not causally connected: {} {}", cb, b);
      return false;
    }

  }
  LOG_TRACE("Adding Blib to emerging Cluster");
  c.blibs_.insert(c.blibs_.end(), b);

  LOG_DEBUG("Leaving TryInsertHit_Emergence()");
  return true;
}


template <class tBlib>
bool TDBScan_Algo<tBlib>::TryInsertHit_Established(
  CausalCluster<tBlib>& c,
  const tBlib& b) {
  LOG_DEBUG("Entering TryInsertHit_Established()");

  auto _result = false;

  int active_connectees = 0;
  auto cb_riter = c.blibs_.crbegin();
  while (cb_riter != c.blibs_.crend()) {
    if (cb_riter->timeTo(b) >= params_.multiplicityTimeWindow) {
      /// we are past the timeframe;
      break;
    }
    active_connectees += CausallyConnected(*cb_riter, b);
    if (active_connectees >= params_.multiplicity) {
      _result = true;
      break;
    }
    ++cb_riter;
  }

  if (_result) {
    LOG_DEBUG("Sufficient Multiplicity in causal overlap; Adding Hit");
    c.blibs_.insert(c.blibs_.end(), b);
    _result = true;
  }

  LOG_DEBUG("Leaving TryInsertHit_Established()");
  return _result;
};



} //namespace tdbscan
