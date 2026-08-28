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
TDBScan_Algo<tBlib>::TDBScan_Algo (
  const TDBScan_ParameterSet& params,
  const Connector_t* connector) :
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
  clusters_.clear();
  active_clusters_.clear();

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
  for (const auto& c : clusters_)
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
  active_clusters_.clear();
  clusters_.clear();

  //process through machinery
  for (const auto& b: blibs){
    log_debug("next Blib");
    NextBlib(b);
  }

  log_debug("Finalize");
  Finalize();

  //prepare output
  BlibSetSequence bss;
  for (const auto& c : clusters_)
    bss.insert(BlibSet(c.hits.begin(), c.hits.end() ));

  log_debug("Leaving Process()");
  return bss;
};


template <class tBlib>
bool TDBScan_Algo<tBlib>::CausallyConnected(const tBlib& b1, const tBlib& b2) const {
  return CausallyConnected(b1, b2, &connector_);
};


template <class tBlib>
void TDBScan_Algo<tBlib>::NextBlib (const tBlib& b) {
  log_debug("Entering NextBlib()");
  //tracking
  bool addedToCluster = false; //keep track of whether h has been added to any cluster

  const auto now = b.GetTime();

  // advance every each cluster in time and try to add the hit to it

  auto ac_iter = active_clusters_.begin();
  while (ac_iter != active_clusters_.end()) {
    //each cluster is advanced in time:
    //removing all too old/expired hits, which cannot make any connections any more;
    //concluded clusters, which do not have any connecting hits left, become 'Inactive' and are put to the garbage
    //if the cluster is still active, try to add the Hit to the cluster

    AdvanceInTime(*ac_iter, now); // check

    if (ac_iter->isConcluded()) {
      if (ac_iter->isEstablished()) {
        clusters_.insert(clusters_.end(), *ac_iter);
      }
      active_clusters_.erase(ac_iter);
    }
    ac_iter++;
  }

  //what remains in the active_cluster_ list are those that can in principle add the hit
  ac_iter = active_clusters_.begin();
  while (ac_iter != active_clusters_.end()) {
    const auto n_overlap = ConnectsTo(*ac_iter, b); //TODO overload function with a fast break after multiplicity is met
    // cluster is established and enough multiplicity -> Grow Cluster
    if (ac_iter->isEstablished()) {
      if (n_overlap >= params_.multiplicity) {
        AddHitToCluster(*ac_iter, b); //<- this needs to internally evaluate (isEstablished)-condition and set it on the cluster
        //else if (n_overlap > 0)
        //  active_clusters_.insert(active_clusters_.end(), ::ConstructSubCluster(*ac_iter, b));
      }
    }
    // cluster is not established yet, but there is some overlap but just not enough
    //else if (!ac_iter.isEstablished()) {

    if (ac_iter->GetCount() == n_overlap) {
      AddHitToCluster(*ac_iter, b); //<- this needs to evaluate (isEstablished)-condition and set it on the cluster
    }
    else {
      //subclusters of those hits overlapping need to be constructed and added to the active_clusters
      active_clusters_.insert(active_clusters_.end(), ConstructSubCluster(*ac_iter, b));
    }
  }

  // at the very end place the hit in a cluster by its own
  active_clusters_.insert(active_clusters_.end(), CausalCluster(b));

  log_debug("Leaving AddHit()");
}




template <class tBlib>
bool TDBScan_Algo<tBlib>::TryInsertHit(
  CausalCluster<tBlib>& c,
  const tBlib& b) {
  log_debug("Entering InsertHit()");

  if (c.status == CausalCluster<tBlib>::CONCLUDED) {
    log_error("cannot add to a concluded cluster")
    return false;
  }

  if (status == CausalCluster<tBlib>::DYING) {
    // or something else

    log_error("cannot add to a dying cluster")
    return false;
  }


  if (status == CausalCluster<tBlib>::EMERGING) {
    int active_connectees = 0;
    for (const auto& cb : c.hits)
      active_connectees += CausallyConnected(cb, b);
    if (active_connectees == c.count())
      // the hit can be added
      return true;
  }

  if (status == CausalCluster<tBlib>::GROWING) {
    int active_connectees = 0;
    for (const auto& cb : c.hits)
      active_connectees += CausallyConnected(cb, b);
    if (active_connectees >= params_.multiplicity) {
      // the hit can be added
      return true;
    }


  c.hits.insert(b);

  if ()

  if (status == Status::DYING) {
    log_error("cannot add to a concluded cluster")
    return;
  }




  log_debug("Leaving InsertHit()");
  return true;
};


void TDBScan_Algo::AddSubEvent(AbsHitSet newSet) {
  log_debug("Entering AddSubEvent()");
  //find any existing subevents which overlap the new one, and merge them into it


}




Time TDBScan_Algo::FinalizedUntil() const {

}

} //namespace tdbscan
