/**
 * \file HiveTrigger.cxx
 *
 * (c) 2012 the IceCube Collaboration
 *
 * $Id: HiveTrigger.cxx 150700 2016-10-12 14:16:27Z mzoll $
 * \version $Revision: 150700 $
 * \date $Date: 2016-10-12 16:16:27 +0200 (ons, 12 okt 2016) $
 * \author Marcel Zoll <marcel.zoll@fysik.su.se>
 */

#include "IceHiveZ/algorithms/HiveTrigger.h"

#include "icetray/I3Units.h"
#include <algorithm>
#include <math.h>
#include <boost/foreach.hpp>

using namespace std;
using namespace HitSorting;

using namespace hivetrigger;

//=============== class HiveTrigger_ParameterSet =================

hivetrigger::HiveTrigger_ParameterSet::HiveTrigger_ParameterSet():
  multiplicity(3),
  multiplicityTimeWindow(1000.*I3Units::ns),
  acceptTimeWindow(NAN),
  rejectTimeWindow(INFINITY),
  connectorBlock(),
  mergeOverlap(1)
{};


//=============== namespace hivetrigger::details =================

inline
bool hivetrigger::detail::CausallyConnected(
  const AbsDAQHit& h1,
  const AbsDAQHit& h2,
  const ConnectorBlockConstPtr& connectorBlock) 
{
  if (h1.GetDAQTicks() > h2.GetDAQTicks())
    return CausallyConnected(h2, h1, connectorBlock); //recursive call to enforce timeorder at this point  
  return connectorBlock->Connected(h1, h2);
}

bool hivetrigger::detail::CausallyOverlaps (
  const AbsDAQHitSet& set1,
  const AbsDAQHitSet& set2,
  const size_t multiplicity,
  const DAQTicks multiplicityTimeWindow)
{
  ///search for identical hits, store them with their hittime, shot them down if passed beyond the timewindow
  typedef std::map<CompactHash, DAQTicks> DOMHitTimes;
  DOMHitTimes common_hitdoms;
  
  AbsDAQHitSet::const_iterator iter1 = set1.begin();
  const AbsDAQHitSet::const_iterator end1 = set1.end();
  AbsDAQHitSet::const_iterator iter2 = set2.begin();
  const AbsDAQHitSet::const_iterator end2 = set2.end();
  while (iter1!=end1 && iter2!=end2) {
    if (*iter1<*iter2)
      ++iter1;
    else if (*iter2<*iter1)
      ++iter2;
    else { //*iter1==*iter2 ; Hits are identical in time and DOM
      const DAQTicks hit_time = iter1->GetDAQTicks();
      //eliminate DOMs where times have run out
      for (DOMHitTimes::iterator it=common_hitdoms.begin(); it!=common_hitdoms.end(); it++) {
        if (it->second<hit_time-multiplicityTimeWindow)
          common_hitdoms.erase(it);
      }
      //add a entry for this DOM
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

hivetrigger::detail::CausalCluster::CausalCluster(
  const HiveTrigger_ParameterSet* p):
  params(p),
  sync_time(std::numeric_limits<DAQTicks>::min()),
  established(false)
{};

inline
DAQTicks hivetrigger::detail::CausalCluster::getEarliestTime() const{
  if (!concluded_hits.empty())
    return(concluded_hits.begin()->GetDAQTicks());
  if (!active_hits.empty())
    return(active_hits.begin()->GetDAQTicks());
  assert(false); //a part of the code, where we should never end up 
  return(std::numeric_limits<DAQTicks>::max());
}

inline
DAQTicks hivetrigger::detail::CausalCluster::getLatestTime() const{
  if (!active_hits.empty())
    return(active_hits.rbegin()->GetDAQTicks());
  if (!concluded_hits.empty())
    return(concluded_hits.rbegin()->GetDAQTicks());
  assert(false); //a part of the code, where we should never end up 
  return(-std::numeric_limits<DAQTicks>::max());
}


bool hivetrigger::detail::CausalCluster::connectsTo(const AbsDAQHit &h) const {
  if (firstHitTimes.count(h.GetDOMIndex())
    && (firstHitTimes.at(h.GetDOMIndex())-h.GetDAQTicks() < NsToTicks(params->acceptTimeWindow)) 
    && (firstHitTimes.at(h.GetDOMIndex())-h.GetDAQTicks() < NsToTicks(params->rejectTimeWindow)))
  {
    return true;
  }
  
  //more elaborate: determine if enough DOMs or all active hits currently in the cluster are connected
  std::set<CompactHash> connectedDOMs; //
  bool allConnected=true;
  for (AbsDAQHitSet::reverse_iterator it=active_hits.rbegin(), end=active_hits.rend(); it!=end; ++it) {
    if (it->GetDOMIndex() == h.GetDOMIndex()) {
      //the DOM of h itself is never to be considered connected
      const DAQTicks dt = h.GetDAQTicks() - it->GetDAQTicks();
      assert(dt>=0);
      if (dt <= NsToTicks(params->acceptTimeWindow)) { // it and h connected
        continue;
      }
      if (dt > NsToTicks(params->rejectTimeWindow) // it rejects h, so not connected
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

void hivetrigger::detail::CausalCluster::insertActiveHit(const AbsDAQHit &h) {
  sync_time = std::max(sync_time, h.GetDAQTicks());
  ++active_doms[h.GetDOMIndex()];
  
  //take care about the first hit-time of any each dom, if this option is enabled
  DOMHitTimes::iterator it= firstHitTimes.find(h.GetDOMIndex());
  if (it==firstHitTimes.end()) //its the first hit on the dom, take its time
    firstHitTimes[h.GetDOMIndex()]=h.GetDAQTicks();
  else if (it->second > h.GetDAQTicks())
    it->second=h.GetDAQTicks();
    
  active_hits.insert(active_hits.end(), h);
  //if the total number of DOMs meets the multiplicity threshold, make note,
  //and also record that this is the last known hit within the cluster contributing
  if (active_doms.size()>=params->multiplicity) {
    established=true;
  }
}


hivetrigger::detail::CausalCluster 
hivetrigger::detail::CausalCluster::getSubCluster(const AbsDAQHit &h) const {
  CausalCluster newSubCluster(params);
  
  if (! firstHitTimes.count(h.GetDOMIndex())) { //FAST short-cut
    //never seen this DOM being hit before; insert them all as long as they are causally conneted
    for (AbsDAQHitSet::iterator it=active_hits.begin(), end=active_hits.end(); it!=end; ++it) {
      if (CausallyConnected(*it, h, params->connectorBlock))
        newSubCluster.insertActiveHit(*it);
    }
  }
  else { //SLOW
    //else: the iteration need to include the check for the accept and rejectTimeWindow
    for (AbsDAQHitSet::iterator it=active_hits.begin(), end=active_hits.end(); it!=end; ++it) {
      if (it->GetDOMIndex() == h.GetDOMIndex()) {
        //check the acceptanceTimeWindow condition
        const DAQTicks dt = h.GetDAQTicks() - it->GetDAQTicks();
        assert(dt >=0); //positive if 'it' earlier than 'h' (the anticipated case)
        if (dt > NsToTicks(params->rejectTimeWindow))
          continue;
        if (dt <= NsToTicks(params->acceptTimeWindow)) {
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
void hivetrigger::detail::CausalCluster::takeConcludedHits (const CausalCluster& c){
  concluded_hits.insert(c.concluded_hits.begin(), c.concluded_hits.end());
}

inline
const AbsDAQHitSet& hivetrigger::detail::CausalCluster::getActiveHits() const {
  return active_hits;
}

inline
const AbsDAQHitSet& hivetrigger::detail::CausalCluster::getConcludedHits() const {
  return concluded_hits;
}

inline
const hivetrigger::detail::CausalCluster::DOMHitTimes& 
hivetrigger::detail::CausalCluster::getFirstHitTimes() const {
  return firstHitTimes;
}

inline
const AbsDAQHit& hivetrigger::detail::CausalCluster::getLatestActiveHit() const{
  return *active_hits.rbegin();
}

bool hivetrigger::detail::CausalCluster::isActive() const {
  if (!active_hits.empty())
    return true;
  else {
    if (params->acceptTimeWindow <= params->multiplicityTimeWindow) {
      //then only active hits can accept more hits
      return false;
    } 
    else {
      //need to look into the firsthit-times if any DOM can still accept a new hit within the acceptanceTimeWindow
      BOOST_FOREACH (const DOMHitTimes::value_type& fht, firstHitTimes) {
        if (fht.second > NsToTicks(sync_time-params->acceptTimeWindow))
          return true;
      }
    }
    return false;
  }
}

inline
bool hivetrigger::detail::CausalCluster::isEstablished() const {
  return established;
}

bool hivetrigger::detail::CausalCluster::isSubsetOf(
  const CausalCluster& c2) const
{
  if (c2.active_hits.size()<this->active_hits.size())
    return(false);
  //use the fact that strict timeorder is enforced in the .active_hits
  AbsDAQHitSet::const_iterator it1=this->active_hits.begin();
  AbsDAQHitSet::const_iterator end1=this->active_hits.end();
  AbsDAQHitSet::const_iterator it2=c2.active_hits.begin();
  AbsDAQHitSet::const_iterator end2=c2.active_hits.end();
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


void hivetrigger::detail::CausalCluster::advanceInTime (
  const DAQTicks ticks) 
{
  while (!active_hits.empty()) {
    const AbsDAQHitSet::const_iterator h=active_hits.begin();
    if (ticks > h->GetDAQTicks()+ NsToTicks(params->multiplicityTimeWindow)) {//the hit is no longer active

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
          BOOST_FOREACH(const AbsDAQHit& hh, active_hits) {
            if (h->GetDOMIndex() == hh.GetDOMIndex())
              firstHitTimes[h->GetDOMIndex()] = hh.GetDAQTicks();
          }
          //NOTE by this shift some inconsitency is introduced of the connections between hits in the cluster
          //however the merging of Clusters in the HiveTrigger will bring this all in sync again
        }
      }
      active_hits.erase(h);
    }
    else
      break;
  }
  //the cluster is now synced to this time
  sync_time=ticks;
}


//===============class HiveTrigger=================================

using namespace hivetrigger::detail;

HiveTrigger::HiveTrigger (const hivetrigger::HiveTrigger_ParameterSet& params):
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
  
  log_info("This is HiveTrigger!");
  log_debug("Leaving Init()");
};

void HiveTrigger::AddHit (const AbsDAQHit& h) {
  log_debug("Entering AddHit()");
  newClusters_.clear();
  bool addedToCluster=false; //keep track of whether h has been added to any cluster

  CausalClusterList::iterator cluster=clusters_.begin();
  while (cluster != clusters_.end()) {
    //each cluster is advanced in time:
    //removing all too old/expired hits, which cannot make any connections any more;
    //concluded clusters, which do not have any connecting hits left, become 'Inactive' and are put to the garbage
    //if the cluster is still active, try to add the Hit to the cluster
    cluster->advanceInTime(h.GetDAQTicks());
    
    if (cluster->isActive()) {
      addedToCluster |= AddHitToCluster(*cluster, h);
      ++cluster;
    }
    else { // if (!cluster->isActive())
      if (cluster->isEstablished()) {
        AbsDAQHitSet subev= cluster->getConcludedHits();
        cluster = clusters_.erase(cluster);
        AddSubEvent(subev);
      }
      else
        cluster = clusters_.erase(cluster);
    }
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


bool HiveTrigger::AddHitToCluster (
  CausalCluster& c,
  const AbsDAQHit& h)
{
  log_debug("Entering AddhitToCluster()");
  if (c.getFirstHitTimes().count(h.GetDOMIndex())
    && (c.getFirstHitTimes().at(h.GetDOMIndex())-h.GetDAQTicks() < NsToTicks(params_.acceptTimeWindow)) 
    && (c.getFirstHitTimes().at(h.GetDOMIndex())-h.GetDAQTicks() < NsToTicks(params_.rejectTimeWindow)))
  {
    c.insertActiveHit(h);
    return true;
  }
  
  
  //more elaborate: determine if enough DOMs or all active hits currently in the cluster are connected
  std::set<CompactHash> connectedDOMs;
  AbsDAQHitSet connectedHits;
  bool allConnected=true;
  
  AbsDAQHitSet::const_reverse_iterator it=c.getActiveHits().rbegin();
  const AbsDAQHitSet::const_reverse_iterator end=c.getActiveHits().rend();
  
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
        const DAQTicks dt = h.GetDAQTicks() - it->GetDAQTicks();
        assert(dt >=0); //h should always be the latest hit
        if (dt <= NsToTicks(params_.acceptTimeWindow)) { // it and h connected
          connectedHits.insert(connectedHits.begin(), *it);
          continue;
        }
        
        if (dt > NsToTicks(params_.rejectTimeWindow)) { // it rejects h, so not connected
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
  BOOST_FOREACH(const AbsDAQHit& connectedHit, connectedHits)
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


void HiveTrigger::AddSubEvent(AbsDAQHitSet& newSet) {
  log_debug("Entering AddSubEvent()");
  
  //find any existing subevents which overlap the new one, and merge them into it
  AbsDAQHitSetList::iterator set =partialSubEvents_.begin();
  while (set != partialSubEvents_.end()) {
    //determine if the overlap sufficent: common hits on 'params.mergeOverlap' DOMs within the time-window
    const bool sufficent_overlap = CausallyOverlaps(newSet,
                                                    *set,
                                                    params_.mergeOverlap,
                                                    NsToTicks(params_.multiplicityTimeWindow));
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
  DAQTicks earliestUpcomingTime = std::numeric_limits<DAQTicks>::max();
  BOOST_FOREACH(const CausalCluster &cluster, clusters_)
    earliestUpcomingTime=std::min(earliestUpcomingTime,cluster.getEarliestTime());

  //any partial subevent whose last hit time is before the earliest time found above
  //cannot be merged again, and so is complete  
  if (earliestUpcomingTime!=std::numeric_limits<DAQTicks>::max()) {
    HiveTrigger::PushEvents(earliestUpcomingTime);
  }
}

void HiveTrigger::PushEvents(const DAQTicks earliestTick) {
  AbsDAQHitSetList::iterator hset=partialSubEvents_.begin();
  while (hset != partialSubEvents_.end()) {
    if (hset->rbegin()->GetDAQTicks() < earliestTick) {
      subEvents_.insert(subEvents_.end(),*hset);
      hset = partialSubEvents_.erase(hset);
    }
    else
      ++hset;
  }
}


void HiveTrigger::AdvanceTime(const DAQTicks ticks) {
  log_debug("Entering AdvanceTime()");

  CausalClusterList::iterator cluster=clusters_.begin();
  while (cluster != clusters_.end()) {
    //each cluster is advanced in time:
    //moving all too old/expired hits, which cannot make any connections any more, out of the active window,
    //concluded clusters, which do not have any active hits left, become 'Inactive' and are put to the garbage
    //or are, in case they are etablished, made into a subevent
    cluster->advanceInTime(ticks);
    
    if (cluster->isActive()) {
      ++cluster;
    }
    else { // if (!cluster->isActive())
      if (cluster->isEstablished()) {
        AbsDAQHitSet subev= cluster->getConcludedHits();
        cluster = clusters_.erase(cluster);
        AddSubEvent(subev);
      }
      else
        cluster = clusters_.erase(cluster);
    }
  }
  log_debug("Leaving AdvanceTime()");  
};

void HiveTrigger::FinalizeSubEvents() {
  log_debug("Entering FinalizeSubEvents()");
  //dump all hits out of the clusters in progress
  
  AdvanceTime(std::numeric_limits<DAQTicks>::max());
  assert(clusters_.size()==0);
  PushEvents(std::numeric_limits<DAQTicks>::max());
  assert(partialSubEvents_.size()==0);
};


AbsDAQHitSetSequence HiveTrigger::PullSubEvents() {
  //return a simple copy of all finisihed subEvents and clear the internal state
  AbsDAQHitSetSequence output = subEvents_;
  subEvents_.clear();
  return output;
};


DAQTicks HiveTrigger::FinalizedUntil() const {
  DAQTicks ticks_frombelow = std::numeric_limits<DAQTicks>::min();
  
  BOOST_FOREACH(const AbsDAQHitSet &sub, subEvents_) {
    //the endtimes of each finished event
    ticks_frombelow = std::max(ticks_frombelow, sub.rbegin()->GetDAQTicks());
  }
  
  DAQTicks ticks_fromabove = std::numeric_limits<DAQTicks>::max();
  BOOST_FOREACH(const AbsDAQHitSet &set, partialSubEvents_) {
    //the endtimes of each finished event
    ticks_frombelow = std::min(ticks_fromabove, set.begin()->GetDAQTicks());
  }
  
  BOOST_FOREACH(const CausalCluster &cluster, clusters_) {
    //the earliest time of all still active clusters
    ticks_fromabove = std::min(ticks_fromabove, cluster.getEarliestTime());
  }
  
  DAQTicks max_time = std::max(ticks_frombelow, ticks_fromabove);
  return (max_time>1 ? max_time-1 : 0); //NOTE 1 because DAQ precision is 1/10ns
}
