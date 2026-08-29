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

#include "tdbscan_algo/dummy_logging.h"

//===========================================
//============== IMPLEMENTATION =============
//===========================================

namespace tdbscan {

//=============== namespace tdbsscan::details =================

//=============== class TDBScan_ParameterSet =================

TDBScan_ParameterSet::TDBScan_ParameterSet():
  multiplicity(3),
  multiplicityTimeWindow( Timediff_t(1000.))
{};


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
  if (params_.acceptTimeWindow<0.0)
    log_fatal("AcceptTimeWindow cannot be negative");
  if (params_.rejectTimeWindow<0.0)
    log_fatal("RejectTimeWindow cannot be negative");
  if (params_.mergeOverlap==0)
    log_warn("RequiredDOMOverlap configured with 0, everything will be merged");

  if (params_.rejectTimeWindow <= params_.acceptTimeWindow)
    log_fatal("RejectTimeWindow needs to be greater than AcceptTimeWindow");

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
  // TODO implement me:
  // something like like: move all established clusters from the active_cluster list into the (nominal) cluster list;
  // observer merging rules
  //sync_time = std::numeric_limits<Time_t>::max();

  assert(active_clusters_.empty());
  return;
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
    log_debug("next Blib");
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
  const auto now = b.getTime();
  // advance every each cluster in time and try to add the hit to it

  // 10. go through all active clusters  and see if blibs have fallen out of the emergence time window and multiplicity cannot be fullfilled; mark them as 'dying'; if there is nothing left mark as 'concluded'
  // 20. go through all emerging clusters and see if blibs have fallen out of the emergence time window, kill them off
  // 30. try to add the hit to existing emerging clusters; if it was added put the emerging clusters on a new_established list;
  // 40. try to add the hit to existing active clusters;
  // 50. traverse the newly established list and try to merge clusters with the active clusters
  // 9. put the hit on a newly created cluster by its own


  // 20. go through all emerging clusters and see if blibs have fallen out of the emergence time window, if so kill the cluster off
  // 21. try to add the hit to remaining clusters; if it was added and cluster establishes (multiplicity met) put the clusters on a new_established list;
  std::list<CausalCluster<tBlib>> _newly_established_clusters;
  auto ec_iter = emerging_clusters_.begin();
  while (ec_iter != emerging_clusters_.end()) {
    const auto _n_active = ec_iter->nHitsWithinTimeWindow(  now, std::numeric_limits<Time_t>::max());// TODO Time type needs proper numeric limits defined
    if (_n_active < ec_iter->count()) {
      ec_iter = emerging_clusters_.erase(ec_iter);
    }

    const auto success = TryInsertHit(*ec_iter, b);
    if (success && ec_iter->count() == params_.multiplicity) {
      _newly_established_clusters.push_back(*ec_iter);
      ec_iter = emerging_clusters_.erase(ec_iter);
      continue;
    }

    ec_iter++;
  }

  // 10. go through all active clusters and see if blibs have fallen out of the emergence time window and multiplicity cannot be fullfilled; mark them as 'dying'; if there is nothing left mark as 'concluded'

  auto ac_iter = active_clusters_.begin();
  while (ac_iter != active_clusters_.end()) {
    const auto _n_active = ac_iter->nHitsWithinTimeWindow(  now, std::numeric_limits<Time_t>::max()); // TODO Time type needs proper numeric limits defined
    if (_n_active == 0) {
      ac_iter->status = CausalCluster<tBlib>::CONCLUDED;
      concluded_clusters_.push_back(*ac_iter);
      ac_iter = active_clusters_.erase(ac_iter);
      continue;
    }
    if (_n_active < params_.multiplicity) {
      ac_iter->status = CausalCluster<tBlib>::DYING; //this is a speedup but algorithmically has no effect
    }
    else {
      const auto success = TryInsertHit(*ec_iter, b);
    }
    ac_iter++;
  }

  // 50. traverse the newly established list and try to merge clusters with the active clusters
  auto nec_iter = _newly_established_clusters.begin();
  ac_iter = active_clusters_.begin();
  while (nec_iter != _newly_established_clusters.end()) {
    while (ac_iter != active_clusters_.end()) {
      if (nec_iter->isSubsetOf(*ac_iter)) {
        nec_iter = _newly_established_clusters.erase(nec_iter);
        ac_iter = active_clusters_.begin();
        continue;
      }
      ac_iter++;
    }
    //established cluster is a genuinely new cluster
    nec_iter->status = CausalCluster<tBlib>::GROWING;
    nec_iter++;
    ac_iter = active_clusters_.begin();
  }

  active_clusters_.insert(active_clusters_.end(), _newly_established_clusters.begin(), _newly_established_clusters.end());
  _newly_established_clusters.clear();

  // 9. put the hit on a newly created cluster by its own
  emerging_clusters_.insert(active_clusters_.end(), CausalCluster(b));

  log_debug("Leaving NextHit()");
}


template <class tBlib>
bool TDBScan_Algo<tBlib>::TryInsertHit(
  CausalCluster<tBlib>& c,
  const tBlib& b) {
  log_debug("Entering InsertHit()");

  switch (c.status) {
    case CausalCluster<tBlib>::CONCLUDED: {
      log_debug("cannot add to a concluded cluster");
      return false;
    }

    case CausalCluster<tBlib>::DYING: {
      // or something else
      log_debug("cannot add to a dying cluster");
      return false;
    }

    case CausalCluster<tBlib>::EMERGING: {
      int active_connectees = 0;
      for (const auto& cb : c.blibs_) {
        if (cb.timeDiff(b) >= params_.emergenceTimeWindow) {
          /// we are past the timeframe;
          return false;
        }
        if (not CausallyConnected(cb, b))
          return false;
      }
      break;
    }

    case CausalCluster<tBlib>::GROWING: {
      int active_connectees = 0;
      for (const auto& cb : c.blibs_) {
        if (cb.timeDiff(b) >= params_.multiplicityTimeWindow)
          /// we are past the timeframe;
            break;
        active_connectees += CausallyConnected(cb, b);
      }
      if (active_connectees >= params_.multiplicity) {
        // the hit can be added
        break;
      }
    }
  }
  c.blibs_.insert(c.blibs_.end(), b);

  log_debug("Leaving TryInsertHit()"); //
  return true;
};



} //namespace tdbscan
