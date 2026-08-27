/**
 * \file HiveCleaning.h
 *
 * (c) 2012 the IceCube Collaboration
 *
 * $Id: HiveCleaning.h 99900 2013-02-26 10:10:43Z mzoll $
 * \version $Revision: 99900 $
 * \date $Date: 2013-02-26 11:10:43 +0100 (Tue, 26 Feb 2013) $
 * \author Marcel Zoll <marcel.zoll@fysik.su.se>
 *
 * The central algorithm to split I3RecoPulseSeriesMaps by arguments of clustering
 */

#ifndef HIVECLEANING_H
#define HIVECLEANING_H

#include <limits>
#include <list>
#include <map>
#include <sstream>

#include "ToolZ/OMKeyHash.h"
#include "ToolZ/HitSorting.h"
#include "IceHiveZ/internals/Connector.h"


/// A set of parameters that steer HiveCleaning
struct HiveCleaning_ParameterSet{
  //parameters
  ///PARAM: A required Multiplicity to the number of surrounding hit DOMs
  unsigned multiplicity;
  /// PARAM: allowed early time residual for the iteration process; an optimization
  double max_tresidual_early;
  /// PARAM: allowed late time residual for the iteration process; an optimization
  double max_tresidual_late;
  /// PARAM: pointer to the ConnectorBlock providing DOM to DOM and Hit connections
  ConnectorBlockPtr connectorBlock;
  
  ///constructor
  HiveCleaning_ParameterSet();
};

///The main cleaning class
class HiveCleaning {
  SET_LOGGER("HiveCleaning");
  
protected://parameters
  //========================
  // Configurable Parameters
  //========================
  /// A parameter-set to run on
  HiveCleaning_ParameterSet params_;

public://methods
  //================
  // Main Interface
  //================
  /// Constructor from a ParameterSet
  HiveCleaning(const HiveCleaning_ParameterSet& params);

  /** @brief ACTION
   * @param hits the hits to process on
   * @return a EventStartStops
   */
  AbsHitSet Clean(const AbsHitSet &hits);
  
  /** @brief ACTION; time order hits first
   * @param hits the hits to process on
   * @return a EventStartStops
   */
  template <class AbsHitContainer>
  AbsHitSet Clean(const AbsHitContainer &hits);
};


//===========================================
//============== IMPLEMENTATION =============
//===========================================

///specialization for already time-ordered hit Containers
template <class AbsHitContainer>
AbsHitSet HiveCleaning::Clean(const AbsHitContainer& inhits) {
  AbsHitSet hs; //timesorted 
  BOOST_FOREACH(const AbsHit& h, inhits)
    hs.insert(h);

  return Clean(hs);
};


#endif //HIVECLEANING_H