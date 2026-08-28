//
// Created by netsu on 28/08/2026.
//



template <class tBlib>
void TDBScan_Algo<tBlib>::NextBlib (const tBlib& h) {
  log_debug("Entering NextBlib()");
  //tracking
  bool addedToCluster = false; //keep track of whether h has been added to any cluster

  const auto now = h.GetTime();

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
    const auto n_overlap = ConnectsTo(*ac_iter, h); //TODO overload function with a fast break after multiplicity is met
    // cluster is established and enough multiplicity -> Grow Cluster
    if (ac_iter->isEstablished()) {
      if (n_overlap >= params_.multiplicity) {
        AddHitToCluster(*ac_iter, h); //<- this needs to internally evaluate (isEstablished)-condition and set it on the cluster
        //else if (n_overlap > 0)
        //  active_clusters_.insert(active_clusters_.end(), ::ConstructSubCluster(*ac_iter, h));
      }
    }
    // cluster is not established yet, but there is some overlap but just not enough
    //else if (!ac_iter.isEstablished()) {

    if (ac_iter->GetCount() == n_overlap) {
      AddHitToCluster(*ac_iter, h); //<- this needs to evaluate (isEstablished)-condition and set it on the cluster
    }
    else {
      //subclusters of those hits overlapping need to be constructed and added to the active_clusters
      active_clusters_.insert(active_clusters_.end(), ConstructSubCluster(*ac_iter, h));
    }
  }

  // at the very end place the hit in a cluster by its own
  active_clusters_.insert(active_clusters_.end(), CausalCluster(h));

  log_debug("Leaving AddHit()");
}


template <class tBlib>
void TDBScan_Algo<tBlib>::AddHitToCluster (
  CausalCluster<tBlib>& c,
  const tBlib& b)
{
  log_debug("Entering AddhitToCluster()");




  if (c.getFirstHitTimes().count(b.GetDOMIndex())
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


void TDBScan_Algo::AddSubEvent(AbsHitSet newSet) {
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


void TDBScan_Algo::FinalizeSubEvents() {
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


Time TDBScan_Algo::FinalizedUntil() const {
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