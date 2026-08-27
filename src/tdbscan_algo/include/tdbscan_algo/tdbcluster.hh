//
// Created by netsu on 08/08/2026.
//


#include "tdbscan_algo/tdbcluster.h"

#include <cstdint>
#include <limits>

namespace tdbscan {
  namespace detail {

template <class tBlib>
bool CausallyConnected(
  const tBlib& h1,
  const tBlib& h2,
  const Connector<tBlib>& connector)
    {
      if (h1.GetTime() > h2.GetTime())
        return CausallyConnected(h2, h1, connector); //recursive call to enforce time-order at this point
      return connector.Connected(h1, h2);
    }


  } // namespace tdbscan::detail


template <class tBlib> 
CausalCluster<tBlib>::CausalCluster() :
  sync_time(-std::numeric_limits<Time_t>::infinity()),
  established(false)
{};

template <class tBlib>
Time_t CausalCluster<tBlib>::getEarliestTime() const {
  if (!hits.empty())
    return(hits.begin()->GetTime());
  return(std::numeric_limits<Time_t>::infinity());
}

template <class tBlib>
Time_t CausalCluster<tBlib>::getLatestTime() const{
  if (!hits.empty())
    return(hits.rbegin()->GetTime());
  return(-std::numeric_limits<Time_t>::infinity());
}

template <class tBlib>
void CausalCluster<tBlib>::insertBlib(const tBlib &h) {
  sync_time = std::max(sync_time, h.GetTime());

  hits.insert(hits.end(), h);
}


template <class tBlib>
void CausalCluster<tBlib>::copyHits (const CausalCluster<tBlib>& c){
  hits.insert(c.hits.begin(), c.hits.end());
}

template <class tBlib>
const typename CausalCluster<tBlib>::BlibSet& CausalCluster<tBlib>::getHits() const {
  return hits;
}


template <class tBlib>
bool CausalCluster<tBlib>::isEstablished() const {
  return established;
}

template <class tBlib>
bool CausalCluster<tBlib>::isConcluded() const {
  return concluded;
}

template <class tBlib>
bool CausalCluster<tBlib>::empty() const {return hits.empty();}

template <class tBlib>
uint64_t CausalCluster<tBlib>::count() const {return hits.count();}



template <class tBlib>
bool CausalCluster<tBlib>::isSubsetOf(
  const CausalCluster& c2) const
{
  if (c2.hits.size()<this->hits.size())
    return(false);
  //use the fact that strict time-order is enforced on the hits
  auto it1=this->hits.begin();
  auto end1=this->hits.end();
  auto it2=c2.hits.begin();
  auto end2=c2.hits.end();
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

};
