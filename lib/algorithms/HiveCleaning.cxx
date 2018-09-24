/**
 * \file HiveCleaning.cxx
 *
 * (c) 2012 the IceCube Collaboration
 *
 * $Id: HiveCleaning.cxx 153493 2017-02-23 17:13:21Z mzoll $
 * \version $Revision: 153493 $
 * \date $Date: 2017-02-23 18:13:21 +0100 (Thu, 23 Feb 2017) $
 * \author Marcel Zoll <marcel.zoll@fysik.su.se>
 */

#include "IceHiveZ/algorithms/HiveCleaning.h"

#include "icetray/I3Units.h"

#include <math.h>
#include <algorithm>
#include <boost/foreach.hpp>

using namespace std;
using namespace HitSorting;


//===============class HiveCleaning_ParameterSet ===================

HiveCleaning_ParameterSet::HiveCleaning_ParameterSet():
  multiplicity(1),
  max_tresidual_early(-INFINITY),
  max_tresidual_late(INFINITY)
{}


//===============class HiveCleaning=================================

HiveCleaning::HiveCleaning(const HiveCleaning_ParameterSet& params) :
  params_(params)
{};


AbsHitSet HiveCleaning::Clean (const AbsHitSet& hits) {
  log_debug("Entering Clean()");
  using namespace HitSorting;

  AbsHitSet outhits;

  if (hits.size()==0) {
    log_warn("The series of hits is empty; Will do nothing");
    return outhits;
  }

  log_debug("Starting Cleaning routine");
  for (AbsHitSet::const_iterator hit_iter = hits.begin(); hit_iter !=hits.end(); ++hit_iter) { //for all hits
    log_trace_stream(" Probing next hit: " << *hit_iter);
    size_t connected_neighbors=0;
    AbsHitSet::const_reverse_iterator past_riter(hit_iter);
    while (past_riter != hits.rend() //
      && (hit_iter->GetTime() - past_riter->GetTime())<=params_.max_tresidual_early)
    { // iterate over all past hits within the time limitation
      if (params_.connectorBlock->Connected(*hit_iter, *past_riter)) { //find connected neighbours
        log_trace_stream("found a past hit to link to : " << *past_riter);
        ++connected_neighbors;
      }
      ++past_riter; //try the next possible neighbour
    }
    
    AbsHitSet::const_iterator future_iter(hit_iter);
    while (future_iter!=hits.end()
      && (future_iter->GetTime() - hit_iter->GetTime())<=params_.max_tresidual_late)
    {
      if (params_.connectorBlock->Connected(*future_iter, *hit_iter)) {
        log_trace_stream("found a future hit to link to : " << *future_iter);
        ++connected_neighbors;
      }
      ++future_iter; //try the next possible neighbour
    }
    
    if (connected_neighbors>=params_.multiplicity) {
      log_debug("found enough connected neighbors");
      outhits.insert(outhits.end(), *hit_iter); //and keep the hit
    }  
  }
  log_debug("Finished Cleaning routine");

  log_debug("Leaving Clean()");
  return outhits;
};

