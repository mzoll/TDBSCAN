//
// Created by netsu on 08/08/2026.
//


#include "include/tdbcluster.h"

#include <limits>

template <class tBlib> 
tdbscan::CausalCluster<tBlib>::CausalCluster() :
  sync_time(-std::numeric_limits<Time_t>::infinity()),
  established(false)
{};

template <class tBlib>
tdbscan::Time_t tdbscan::CausalCluster<tBlib>::getEarliestTime() const {
  if (!hits.empty())
    return(hits.begin()->GetTime());
  return(std::numeric_limits<tdbscan::Time_t>::infinity());
}

template <class tBlib>
tdbscan::Time_t tdbscan::CausalCluster<tBlib>::getLatestTime() const{
  if (!hits.empty())
    return(hits.rbegin()->GetTime());
  return(-std::numeric_limits<tdbscan::Time_t>::infinity());
}

template <class tBlib>
void tdbscan::CausalCluster<tBlib>::insertBlib(const tBlib &h) {
  sync_time = std::max(sync_time, h.GetTime());

  hits.insert(hits.end(), h);
}


template <class tBlib>
void tdbscan::CausalCluster<tBlib>::copyHits (const CausalCluster<tBlib>& c){
  hits.insert(c.hits.begin(), c.hits.end());
}

inline
const AbsHitSet& tdbscan::CausalCluster<tBlib>::getActiveHits() const {
  return active_hits;
}

inline
const AbsHitSet& tdbscan::CausalCluster<tBlib>::getConcludedHits() const {
  return concluded_hits;
}

inline
const tdbscan::CausalCluster<tBlib>::DOMHitTimes& 
tdbscan::CausalCluster<tBlib>::getFirstHitTimes() const {
  return firstHitTimes;
}

inline
const AbsHit& tdbscan::CausalCluster<tBlib>::getLatestActiveHit() const{
  return *active_hits.rbegin();
}

bool tdbscan::CausalCluster<tBlib>::isActive() const {
  if (!active_hits.empty())
    return true;

  if (params->acceptTimeWindow <= params->multiplicityTimeWindow) { // TODO watch out; this is outside info
    //then only active hits can accept more hits
    return false;
  }
  //need to look into the time of every first hit on any each DOM if further hits can be accepted
  const auto latest_accept_time = sync_time - params->acceptTimeWindow; // TODO watch out; this is outside info
  BOOST_FOREACH (const DOMHitTimes::value_type& fht, firstHitTimes) {
    if (fht.second > latest_accept_time)
      return true;
  }
  return false;
}

inline
bool tdbscan::CausalCluster<tBlib>::isEstablished() const {
  return established;
}

bool tdbscan::CausalCluster<tBlib>::isSubsetOf(
  const CausalCluster& c2) const
{
  if (c2.active_hits.size()<this->active_hits.size())
    return(false);
  //use the fact that strict timeorder is enforced on the active_hits
  AbsHitSet::const_iterator it1=this->active_hits.begin();
  AbsHitSet::const_iterator end1=this->active_hits.end();
  AbsHitSet::const_iterator it2=c2.active_hits.begin();
  AbsHitSet::const_iterator end2=c2.active_hits.end();
  for (; it1!=end1 && it2!=end2; ++it1, ++it2) {
    //if the two current items don't match scan though the (potential) superset looking for a match
    while (it2!=end2 && *it2<*it1)
      ++it2;
    //three possible cases arise:
    //if the items are now equal, c1 still appears to be a subset
    //if the item in c2 is greater than the one in c1, or we've gone off of the end of c2
    // there is no match for this item, so c1 is not a subset
    if (it2==end2 || *it1<*it2)
      return(false);
  }
  //if all of the items matched until we ran off of the end of c2,
  //but there are still items left in c1, c1 is not a subset
  return(!(it1!=end1 && it2==end2));
}


void tdbscan::CausalCluster<tBlib>::advanceInTime (
  const Time time) 
{
  while (!active_hits.empty()) {
    const auto ah_iter = active_hits.begin();
    if (time > ah_iter->GetTime() + params->multiplicityTimeWindow) {// TODO watch out; this is outside info
      //the hit is no longer active, thus
      //decrement the number of hits on the DOM where h occurred
      if ((--active_doms[ah_iter->GetDOMIndex()])<=0) //NOTE TODO do we need to bother with this after the cluster is established, and this is probably not checked anymore?
        active_doms.erase(ah_iter->GetDOMIndex());

      //if the multiplicity threshold was met include h in the finished cluster
      if (established) {
        //insert the hit
        concluded_hits.insert(concluded_hits.end(),*ah_iter);
      }
      else { //hit is about to be discarded
        //sync up the firsthit-time map
        if (!active_doms.count(ah_iter->GetDOMIndex()))
          firstHitTimes.erase(ah_iter->GetDOMIndex());
        else {
          //check for the next hit on the same DOM which is still active and take its time instead
          BOOST_FOREACH(const AbsHit& hh, active_hits) {
            if (ah_iter->GetDOMIndex() == hh.GetDOMIndex())
              firstHitTimes[ah_iter->GetDOMIndex()] = hh.GetTime();
          }
          //NOTE by this shift some inconsistency is introduced of the connections between hits in the cluster
          //however the merging of Clusters in the HiveSplitter will bring this all in sync again
        }
      }
      active_hits.erase(ah_iter);
    }
    else
      break;
  }

  //the cluster is now synced to this time
  sync_time=time;
}
