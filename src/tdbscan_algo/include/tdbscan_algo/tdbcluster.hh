//
// Created by netsu on 08/08/2026.
//


#pragma once

#include "tdbscan_algo/tdbcluster.h"

#include <cstdint>
#include <limits>

namespace tdbscan {
template <class tBlib>
CausalCluster<tBlib>::CausalCluster()
{};

template <class tBlib>
CausalCluster<tBlib>::CausalCluster(const tBlib &b)
{blibs_.insert(b);};

template <class tBlib>
CausalCluster<tBlib>::CausalCluster(const std::set<tBlib> &bset)
  {blibs_.insert(bset.cbegin(), bset.cend());};

template <class tBlib>
CausalCluster<tBlib>::CausalCluster(const CausalCluster& cc)
{blibs_.insert(cc.blibs_.cbegin(), cc.blibs_.cend());};

template <class tBlib>
Time_t CausalCluster<tBlib>::getEarliestTime() const {
  if (!blibs_.empty())
    return(blibs_.begin()->GetTime());
  return(std::numeric_limits<Time_t>::infinity());
}

template <class tBlib>
Time_t CausalCluster<tBlib>::getLatestTime() const{
  if (!blibs_.empty())
    return(blibs_.rbegin()->GetTime());
  return(-std::numeric_limits<Time_t>::infinity());
}


template <class tBlib>
uint64_t CausalCluster<tBlib>::nHitsWithinTimeWindow(
  const Time_t earliest, const Time_t latest) const {
  uint64_t _count = 0;
  for (const auto& b : blibs_) {
    if (earliest <= b.getTime() && b.getTime() >= latest)
      _count++;
  }
  return _count;
}


template <class tBlib>
void CausalCluster<tBlib>::insertBlib(const tBlib &h) {
  blibs_.insert(blibs_.end(), h);
}


template <class tBlib>
void CausalCluster<tBlib>::copyHits (const CausalCluster<tBlib>& c){
  blibs_.insert(c.blibs_.begin(), c.blibs_.end());
}

template <class tBlib>
const typename CausalCluster<tBlib>::BlibSet& CausalCluster<tBlib>::getHits() const {
  return blibs_;
}

template <class tBlib>
bool CausalCluster<tBlib>::isSubsetOf(
  const CausalCluster& c2) const
{
  if (c2.blibs_.size()<this->blibs_.size())
    return(false);
  //use the fact that strict time-order is enforced on the blibs_
  auto it1=this->blibs_.begin();
  auto end1=this->blibs_.end();
  auto it2=c2.blibs_.begin();
  auto end2=c2.blibs_.end();
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

template <class tBlib>
bool CausalCluster<tBlib>::isSupersetOf(
  const CausalCluster& c2) const
{
  return c2.isSubsetOf(*this);
}

template <class tBlib>
bool CausalCluster<tBlib>::isConcruent(
  const CausalCluster& c2) const
{
  if (c2.blibs_.size() != this->blibs_.size())
    return false;

  // simultaneous step through all members and check they are equal
  auto it1=this->blibs_.begin();
  auto end1=this->blibs_.end();
  auto it2=c2.blibs_.begin();
  auto end2=c2.blibs_.end();
  while (it1 != end1 && it2 != end2) {
    if (it1 != it2)
      return false;
    ++it1;
    ++it2;
  }
  return true;
}

template <class tBlib>
bool CausalCluster<tBlib>::empty() const {return blibs_.empty();}

template <class tBlib>
uint64_t CausalCluster<tBlib>::count() const {return blibs_.size();}

};
