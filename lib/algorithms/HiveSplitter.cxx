/**
 * \file HiveSplitter.cxx
 *
 * (c) 2012 the IceCube Collaboration
 *
 * $Id: HiveSplitter.cxx 153493 2017-02-23 17:13:21Z mzoll $
 * \version $Revision: 153493 $
 * \date $Date: 2017-02-23 18:13:21 +0100 (Thu, 23 Feb 2017) $
 * \author Marcel Zoll <marcel.zoll@fysik.su.se>
 */

#include "IceHiveZ/algorithms/HiveSplitter.h"

#include "icetray/I3Units.h"
#include <algorithm>
#include <math.h>
#include <boost/foreach.hpp>

using namespace std;
using namespace HitSorting;

using namespace hivesplitter;

//=============== class HiveSplitter_ParameterSet =================

hivesplitter::HiveSplitter_ParameterSet::HiveSplitter_ParameterSet():
  multiplicity(3),
  multiplicityTimeWindow(1000.*I3Units::ns),
  acceptTimeWindow(NAN),
  rejectTimeWindow(INFINITY),
  connectorBlock(),
  mergeOverlap(1)
{};


//=============== namespace hivesplitter::details =================

inline
bool hivesplitter::detail::CausallyConnected(
  const AbsHit& h1,
  const AbsHit& h2,
  const ConnectorBlockConstPtr& connectorBlock) 
{
  if (h1.GetTime() > h2.GetTime())
    return CausallyConnected(h2, h1, connectorBlock); //recursive call to enforce timeorder at this point  
  return connectorBlock->Connected(h1, h2);
}

bool hivesplitter::detail::CausallyOverlaps (
  const AbsHitSet& set1,
  const AbsHitSet& set2,
  const size_t multiplicity,
  const Time multiplicityTimeWindow)
{
  ///search for identical hits, store them with their hittime, shot them down if passed beyond the timewindow
  typedef std::map<CompactHash, Time> DOMHitTimes;
  DOMHitTimes common_hitdoms;
  
  AbsHitSet::const_iterator iter1 = set1.begin();
  const AbsHitSet::const_iterator end1 = set1.end();
  AbsHitSet::const_iterator iter2 = set2.begin();
  const AbsHitSet::const_iterator end2 = set2.end();
  while (iter1!=end1 && iter2!=end2) {
    if (*iter1<*iter2)
      ++iter1;
    else if (*iter2<*iter1)
      ++iter2;
    else { //*iter1==*iter2 ; Hits are identical in time and DOM
      const Time hit_time = iter1->GetTime();
      //eliminate DOMs where times have run out
      for (DOMHitTimes::iterator it=common_hitdoms.begin(); it!=common_hitdoms.end(); it++) {
        if (it->second<hit_time-multiplicityTimeWindow)
          common_hitdoms.erase(it);
      }
      common_hitdoms[iter1->GetDOMIndex()] = hit_time;
      if (common_hitdoms.size()>=multiplicity) {
        //found enough required overlap within the time-window
        return true;
      }
      ++iter1;
      ++iter2;
    }
  }
  return false;
};

hivesplitter::detail::CausalCluster::CausalCluster(
  const HiveSplitter_ParameterSet* p):
  params(p),
  sync_time(-std::numeric_limits<Time>::infinity()),
  established(false)
{};

inline
Time hivesplitter::detail::CausalCluster::getEarliestTime() const{
  if (!concluded_hits.empty())
    return(concluded_hits.begin()->GetTime());
  if (!active_hits.empty())
    return(active_hits.begin()->GetTime());
  assert(false); //a part of the code, where we should never end up 
  return(std::numeric_limits<Time>::infinity());
}

inline
Time hivesplitter::detail::CausalCluster::getLatestTime() const{
  if (!active_hits.empty())
    return(active_hits.rbegin()->GetTime());
  if (!concluded_hits.empty())
    return(concluded_hits.rbegin()->GetTime());
  assert(false); //a part of the code, where we should never end up 
  return(-std::numeric_limits<Time>::infinity());
}


bool hivesplitter::detail::CausalCluster::connectsTo(const AbsHit &h) const {
  if (firstHitTimes.count(h.GetDOMIndex())
    && (firstHitTimes.at(h.GetDOMIndex())-h.GetTime() < params->acceptTimeWindow) 
    && (firstHitTimes.at(h.GetDOMIndex())-h.GetTime() < params->rejectTimeWindow))
  {
    return true;
  }
  
  //more elaborate: determine if enough DOMs or all active hits currently in the cluster are connected
  std::set<CompactHash> connectedDOMs; //
  bool allConnected=true;
  for (AbsHitSet::reverse_iterator it=active_hits.rbegin(), end=active_hits.rend(); it!=end; ++it) {
    if (it->GetDOMIndex() == h.GetDOMIndex()) {
      //the DOM of h itself is never to be considered connected
      const Time dt = h.GetTime() - it->GetTime();
      if (dt <= params->acceptTimeWindow) { // it and h connected
        continue;
      }
      if (dt > params->rejectTimeWindow // it rejects h, so not connected
        || ! CausallyConnected(*it, h, params->connectorBlock)) // it cannnot even connect to h
      {
        allConnected = false;
      }
    }
    
    if (CausallyConnected(*it, h, params->connectorBlock)) { //not on the same DOM
      //try if h can connect to hits on other DOMs
      connectedDOMs.insert(it->GetDOMIndex()); // add to the number of connected DOMs
      if (connectedDOMs.size() >= params->multiplicity-1) // found enough connections
        return true;
    }
    else //none of the possible connections of h to it worked out
      allConnected=false;
  }
  return allConnected;
};

void hivesplitter::detail::CausalCluster::insertActiveHit(const AbsHit &h) {
  sync_time = std::max(sync_time, h.GetTime());
  ++active_doms[h.GetDOMIndex()];
  
  //take care about the first hit-time of any each dom, if this option is enabled
  DOMHitTimes::iterator it= firstHitTimes.find(h.GetDOMIndex());
  if (it==firstHitTimes.end()) //its the first hit on the dom, take its time
    firstHitTimes[h.GetDOMIndex()]=h.GetTime();
  else if (it->second > h.GetTime())
    it->second=h.GetTime();
    
  active_hits.insert(active_hits.end(), h);
  //if the total number of DOMs meets the multiplicity threshold, make note,
  //and also record that this is the last known hit within the cluster contributing
  if (active_doms.size()>=params->multiplicity) {
    established=true;
  }
}


hivesplitter::detail::CausalCluster 
hivesplitter::detail::CausalCluster::getSubCluster(const AbsHit &h) const {
  CausalCluster newSubCluster(params);
  
  if (! firstHitTimes.count(h.GetDOMIndex())) { //FAST short-cut
    //never seen this DOM being hit before; insert them all as long as they are causally conneted
    for (AbsHitSet::iterator it=active_hits.begin(), end=active_hits.end(); it!=end; ++it) {
      if (CausallyConnected(*it, h, params->connectorBlock))
        newSubCluster.insertActiveHit(*it);
    }
  }
  else { //SLOW
    //else: the iteration need to include the check for the accept and rejectTimeWindow
    for (AbsHitSet::iterator it=active_hits.begin(), end=active_hits.end(); it!=end; ++it) {
      if (it->GetDOMIndex() == h.GetDOMIndex()) {
        //check the acceptanceTimeWindow condition
        const Time dt = h.GetTime() - it->GetTime(); //positive if 'it' earlier than 'h' (the anticipated case)
        if (dt > params->rejectTimeWindow)
          continue;
        if (dt >=0 && dt <= params->acceptTimeWindow) {
          newSubCluster.insertActiveHit(*it);
          continue;
        }
      }
      if (CausallyConnected(*it, h, params->connectorBlock))
        newSubCluster.insertActiveHit(*it);
    }
  }
  return newSubCluster;
}

inline
void hivesplitter::detail::CausalCluster::takeConcludedHits (const CausalCluster& c){
  concluded_hits.insert(c.concluded_hits.begin(), c.concluded_hits.end());
}

inline
const AbsHitSet& hivesplitter::detail::CausalCluster::getActiveHits() const {
  return active_hits;
}

inline
const AbsHitSet& hivesplitter::detail::CausalCluster::getConcludedHits() const {
  return concluded_hits;
}

inline
const hivesplitter::detail::CausalCluster::DOMHitTimes& 
hivesplitter::detail::CausalCluster::getFirstHitTimes() const {
  return firstHitTimes;
}

inline
const AbsHit& hivesplitter::detail::CausalCluster::getLatestActiveHit() const{
  return *active_hits.rbegin();
}

bool hivesplitter::detail::CausalCluster::isActive() const {
  if (!active_hits.empty())
    return true;
  else {
    if (params->acceptTimeWindow<=params->multiplicityTimeWindow) {
      //then only active hits can accept more hits
      return false;
    } 
    else {
      //need to look into the firsthit-times if any DOM can still accept a new hit within the acceptanceTimeWindow
      BOOST_FOREACH (const DOMHitTimes::value_type& fht, firstHitTimes) {
        if (fht.second > sync_time-params->acceptTimeWindow)
          return true;
      }
    }
    return false;
  }
}

inline
bool hivesplitter::detail::CausalCluster::isEstablished() const {
  return established;
}

bool hivesplitter::detail::CausalCluster::isSubsetOf(
  const CausalCluster& c2) const
{
  if (c2.active_hits.size()<this->active_hits.size())
    return(false);
  //use the fact that strict timeorder is enforced in the .active_hits
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


void hivesplitter::detail::CausalCluster::advanceInTime (
  const Time time) 
{
  while (!active_hits.empty()) {
    const AbsHitSet::const_iterator h=active_hits.begin();
    if (time > h->GetTime()+ params->multiplicityTimeWindow) {//the hit is no longer active

      //decrement the number of hits on the DOM where h occurred
      if ((--active_doms[h->GetDOMIndex()])<=0) //NOTE TODO do we need to bother with this after the cluster is established, and this is probably not checked anymore?
        active_doms.erase(h->GetDOMIndex());

      //if the mutiplicity threshold was met include h in the finished cluster
      if (established) {
        //insert the hit
        concluded_hits.insert(concluded_hits.end(),*h);
      }
      else { //hit is about to be discarded
        //sync up the firsthit-time map
        if (!active_doms.count(h->GetDOMIndex()))
          firstHitTimes.erase(h->GetDOMIndex());
        else {
          //check for the next hit on the same DOM which is still active and take its time instead
          BOOST_FOREACH(const AbsHit& hh, active_hits) {
            if (h->GetDOMIndex() == hh.GetDOMIndex())
              firstHitTimes[h->GetDOMIndex()] = hh.GetTime();
          }
          //NOTE by this shift some inconsitency is introduced of the connections between hits in the cluster
          //however the merging of Clusters in the HiveSplitter will bring this all in sync again
        }
      }
      active_hits.erase(h);
    }
    else
      break;
  }
  //the cluster is now synced to this time
  sync_time=time;
}


//===============class HiveSplitter=================================

using namespace hivesplitter::detail;

HiveSplitter::HiveSplitter (const hivesplitter::HiveSplitter_ParameterSet& params):
  params_(params)
{
  if (params_.multiplicity<=0)
    log_fatal("Multiplicity should be greater than zero");
  if (params_.multiplicityTimeWindow<=0.0)
    log_fatal("TimeWindow should be greater than zero");
  if (params_.acceptTimeWindow<0.0)
    log_fatal("AcceptTimeWindow cannot be negative");
  if (params_.rejectTimeWindow<0.0)
    log_fatal("RejectTimeWindow cannot be negative");
  
  if (params_.rejectTimeWindow <= params_.acceptTimeWindow)
    log_fatal("RejectTimeWindow needs to be greater than AcceptTimeWindow");
  
  if (! params_.connectorBlock)
    log_error("No ConnectionBlock defined!");
  //TODO check integrety of connectorBlock
  
  if (params_.mergeOverlap==0)
    log_warn("RequiredDOMOverlap configured with 0, everything will be merged");
  
  log_info("This is HiveSplitter!");
  log_debug("Leaving Init()");
};

//specialize for AbsHitSet, which is already time-ordered
template <>
AbsHitSetSequence HiveSplitter::Split<AbsHitSet> (const AbsHitSet& inhits) {
  log_debug("Entering Split()");
  clusters_.clear();
  newClusters_.clear();
  partialSubEvents_.clear();
  subEvents_.clear();

  //process through machinery
  BOOST_FOREACH(const AbsHit& h, inhits) {
    log_debug("next Hit");
    AddHit(h);
  }
    
  log_debug("Finalize");
  FinalizeSubEvents();

  log_debug("Leaving Split()");
  return subEvents_;
};

void HiveSplitter::AddHit (const AbsHit& h) {
  log_debug("Entering AddHit()");
  newClusters_.clear();
  bool addedToCluster=false; //keep track of whether h has been added to any cluster

  CausalClusterList::iterator cluster=clusters_.begin();
  while (cluster != clusters_.end()) {
    //each cluster is advanced in time:
    //removing all too old/expired hits, which cannot make any connections any more;
    //concluded clusters, which do not have any connecting hits left, become 'Inactive' and are put to the garbage
    //if the cluster is still active, try to add the Hit to the cluster
    cluster->advanceInTime(h.GetTime());
    
    if (cluster->isEstablished() && !cluster->isActive()) {
      AddSubEvent(cluster->getConcludedHits());
    }
    if (cluster->isActive()) {
      addedToCluster |= AddHitToCluster(*cluster, h);
      ++cluster;
    }
    else //concluded clusters are killed off
      cluster = clusters_.erase(cluster);    
  }

  //Move all newly generated clusters into the main cluster list,
  //eliminating clusters which are subsets of other clusters
  for (CausalClusterList::iterator newCluster=newClusters_.begin(), nend=newClusters_.end(); newCluster!=nend; newCluster++) {
    bool add=true;

    cluster=clusters_.begin();
    while (cluster != clusters_.end()) {
      //check whether the new cluster is a subset of the old cluster
      //if the old cluster does not contain h, it cannot be a superset of the new cluster which does,
      //and if the old cluster contains h, it will be the last hit in that cluster
      if (cluster->getLatestActiveHit()==h){
        if (newCluster->isSubsetOf(*cluster)) {
          add=false;
          break;
        }
        ++cluster;
      }
      //otherwise, the new cluster may still be a superset of the old cluster
      else if (cluster->isSubsetOf(*newCluster)) {
        //if replacing, make sure not to lose any hits already shifted to the old cluster's concluded_hits list
        newCluster->takeConcludedHits(*cluster);
        cluster = clusters_.erase(cluster);
      }
      else
        ++cluster;
    }
    if (add)
      clusters_.push_back(*newCluster);
  }
  newClusters_.clear();

  //if h was not added to any cluster, put it in a cluster by itself
  if (!addedToCluster) {
    clusters_.push_back(CausalCluster(&params_));
    clusters_.back().insertActiveHit(h);
  }
  log_debug("Leaving AddHit()");
}


bool HiveSplitter::AddHitToCluster (
  CausalCluster& c,
  const AbsHit& h)
{
  log_debug("Entering AddhitToCluster()");
  if (c.getFirstHitTimes().count(h.GetDOMIndex())
    && (c.getFirstHitTimes().at(h.GetDOMIndex())-h.GetTime() < params_.acceptTimeWindow) 
    && (c.getFirstHitTimes().at(h.GetDOMIndex())-h.GetTime() < params_.rejectTimeWindow))
  {
    c.insertActiveHit(h);
    return true;
  }
  
  
  //more elaborate: determine if enough DOMs or all active hits currently in the cluster are connected
  std::set<CompactHash> connectedDOMs;
  AbsHitSet connectedHits;
  bool allConnected=true;
  
  AbsHitSet::const_reverse_iterator it=c.getActiveHits().rbegin();
  const AbsHitSet::const_reverse_iterator end=c.getActiveHits().rend();
  
  if (! c.getFirstHitTimes().count(h.GetDOMIndex())) { //FAST
    //never seen the DOM of h being hit before; check just causallyConnected
    for (; it!=end; ++it) {
      if (CausallyConnected(*it, h, params_.connectorBlock)) {
        connectedDOMs.insert(it->GetDOMIndex());
        connectedHits.insert(connectedHits.begin(), *it);
        //exit condition
        if (connectedDOMs.size() >= params_.multiplicity-1) {// found enough connections
          c.insertActiveHit(h);
          return true;
        }
      }
      else {
        allConnected=false;
      }
    }
  }
  else {//SLOW
    //need to check the conditions of accept/reject on same DOM
    for (; it!=end; ++it) {
      if (it->GetDOMIndex() == h.GetDOMIndex()) {
        //the DOM of h itself is never to be considered connected
        const Time dt = h.GetTime() - it->GetTime();
        //assert(dt >=0); //h should always be the latest hit
        if (dt <= params_.acceptTimeWindow) { // it and h connected
          connectedHits.insert(connectedHits.begin(), *it);
          continue;
        }
        
        if (dt > params_.rejectTimeWindow) { // it rejects h, so not connected
          allConnected = false;
          continue;
        }
        
        if ( CausallyConnected(*it, h, params_.connectorBlock))
          connectedHits.insert(connectedHits.begin(), *it);
        else
          allConnected=false;
      }
      else {
        //not on the same DOM
        if (CausallyConnected(*it, h, params_.connectorBlock)) { 
          //try if h can connect to hits on other DOMs
          connectedDOMs.insert(it->GetDOMIndex()); // add to the number of connected DOMs
          connectedHits.insert(connectedHits.begin(), *it);
          //exit condition
          if (connectedDOMs.size() >= params_.multiplicity-1) {// found enough connections
            c.insertActiveHit(h);
            return true;
          }
        }
        else //none of the possible connections of h to it worked out
          allConnected=false;
      }
    }    
  }
  
  if (allConnected) {
    //when all hits, when all hits which are in the cluster are connecting, thats also OKay
    c.insertActiveHit(h);
    return true;
  }
  
  if (connectedHits.empty()) {
    //no overlap at all
    return false;    
  }
  
  CausalCluster newSubCluster(&params_);
  BOOST_FOREACH(const AbsHit& connectedHit, connectedHits)
    newSubCluster.insertActiveHit(connectedHit);
  newSubCluster.insertActiveHit(h); //insert the hit itself now
  
  bool keep=true;
  CausalClusterList::iterator iter=newClusters_.begin();
  while (iter != newClusters_.end()) {
    if (iter->isSubsetOf(newSubCluster))
      iter = newClusters_.erase(iter); //remove a redundant, existing cluster
    else if (newSubCluster.isSubsetOf(*iter)) {
      keep=false; //this cluster is redundant, so abort adding it
      break;
    }
    else
      ++iter;
  }
  if (keep) //finally, actually add the new cluster, as long as it isn't redundant
    newClusters_.push_back(newSubCluster);
  
  return true;
};


void HiveSplitter::AddSubEvent(AbsHitSet newSet) {
  log_debug("Entering AddSubEvent()");
  //find any existing subevents which overlap the new one, and merge them into it
  
  AbsHitSetList::iterator set =partialSubEvents_.begin();
  while (set != partialSubEvents_.end()) {
    //determine if the overlap sufficent: common hits on 'params.mergeOverlap' DOMs within the time-window
    const bool sufficent_overlap = CausallyOverlaps(newSet, *set, params_.mergeOverlap, params_.multiplicityTimeWindow);
    if (sufficent_overlap) {
      newSet.insert(set->begin(),set->end());
      set = partialSubEvents_.erase(set);
    }
    else
      ++set;
  }
  
  partialSubEvents_.push_back(newSet);
  newSet.clear();

  //find the earliest time of all hits currently percolating through the clusters
  Time earliestUpcomingTime = std::numeric_limits<Time>::infinity();
  BOOST_FOREACH(const CausalCluster &cluster, clusters_)
    earliestUpcomingTime=std::min(earliestUpcomingTime,cluster.getEarliestTime());

  //any partial subevent whose last hit time is before the earliest time found above
  //cannot be merged again, and so is complete
  if (earliestUpcomingTime!=std::numeric_limits<Time>::infinity()) {
    AbsHitSetList::iterator hset=partialSubEvents_.begin();
    while (hset != partialSubEvents_.end()) {
      if (hset->rbegin()->GetTime() < earliestUpcomingTime) {
        //copy the contents of this subevent to a new subevent with ordering suitable
        //for retrieval of the actual hits and file it under the time of its first hit
        subEvents_.insert(subEvents_.end(),*hset);
        hset = partialSubEvents_.erase(hset);
      }
      else
        ++hset;
    }
  }
}


void HiveSplitter::FinalizeSubEvents() {
  log_debug("Entering FinalizeSubEvents()");
  //dump all hits out of the clusters in progress
  
  CausalClusterList::iterator cluster = clusters_.begin();
  while (cluster!=clusters_.end()) {
    cluster->advanceInTime(std::numeric_limits<Time>::infinity());
    if (cluster->isEstablished())
      AddSubEvent(cluster->getConcludedHits());
    cluster = clusters_.erase(cluster);
  }
    
  //clusters_.clear(); //should already be empty
  //collect all leftover subevents
  BOOST_FOREACH(AbsHitSet &set, partialSubEvents_)
    subEvents_.insert(subEvents_.end(),set);
  partialSubEvents_.clear();
};


Time HiveSplitter::FinalizedUntil() const {
  Time time_frombelow = -std::numeric_limits<Time>::infinity();
  
  BOOST_FOREACH(const AbsHitSet &sub, subEvents_) {
    //the endtimes of each finished event
    time_frombelow = std::max(time_frombelow, sub.rbegin()->GetTime());
  }
  
  Time time_fromabove = std::numeric_limits<Time>::infinity();
  BOOST_FOREACH(const AbsHitSet &set, partialSubEvents_) {
    //the endtimes of each finished event
    time_frombelow = std::min(time_fromabove, set.begin()->GetTime());
  }
  
  BOOST_FOREACH(const CausalCluster &cluster, clusters_) {
    //the earliest time of all still active clusters
    time_fromabove = std::min(time_fromabove, cluster.getEarliestTime());
  }
  
  Time max_time = std::max(time_frombelow, time_fromabove);
  return max_time-0.1; //NOTE 0.1 because DAQ precision is 1/10ns
}
