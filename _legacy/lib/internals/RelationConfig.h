/**
 * \file HiveRelationConfig.h
 *
 * (c) 2012 the IceCube Collaboration
 *
 * $Id: HiveSplitter.h 99900 2013-02-26 10:10:43Z mzoll $
 * \version $Revision: 99900 $
 * \date $Date: 2013-02-26 11:10:43 +0100 (Tue, 26 Feb 2013) $
 * \author Marcel Zoll <marcel.zoll@fysik.su.se>
 *
 * A confiuration of ringsettings that can be translated into Connectors
 */

#ifndef RELATIONCONFIG_H
#define RELATIONCONFIG_H

#include "IceHiveZ/__SERIALIZATION.h"

static const unsigned relationconfig_version_ = 0;

#include <limits>
#include <list>
#include <map>

#include "dataclasses/I3Constants.h"
#include "ToolZ/GCDinfo.h"
#include "ToolZ/OMKeyHash.h"
#include "ToolZ/Hitclasses.h"
#include "ToolZ/HashedGeometry.h"
#include "IceHiveZ/internals/Relation.h"
#include "IceHiveZ/internals/Hive.h"

#include <boost/make_shared.hpp>
#include <boost/foreach.hpp>


///namespace in order to define parameters for IceHive-services

/**
 * Abstract BasePointer class to derive from
 */
class RelationConfig {
public:
  /// Build a Relation from the information in the derived class
  virtual 
  RelationPtr BuildRelation (
    const HashedGeometryConstPtr& hashedGeo) const=0;
};

typedef boost::shared_ptr<RelationConfig> RelationConfigPtr;
typedef boost::shared_ptr<const RelationConfig> RelationConfigConstPtr;


//===================== CLASS SimpleRelationConfig ==========

/**
 * Encapusules a predication call to construct the Relation
 */
class SimpleRelationConfig : public RelationConfig {
public:
  ///set this object on construction
  boost::function<bool (const OMKey&, const OMKey&)> callobj_;
public: //ctor
  /// constructor
  SimpleRelationConfig(
    boost::function<bool (const OMKey&, const OMKey&)> callobj);
public:
  RelationPtr BuildRelation (const HashedGeometryConstPtr& hashedGeo) const;
};

typedef boost::shared_ptr<SimpleRelationConfig> SimpleRelationConfigPtr;
typedef boost::shared_ptr<const SimpleRelationConfig> SimpleRelationConfigConstPtr;

//===================== CLASS HiveRelationConfiguration =========================

/** Construct a Relation by the Ring and Hive configurations of IceHive
 */
class HiveRelationConfig : public RelationConfig
{
  SET_LOGGER("HiveRelationConfig");
public: //utility classes  
  //===================== STRUCT LimitPairs ==================
  /// pair of two numbers which signal the relative position (in meters) up/down a string 
  /// relative another DOM 
  struct LimitPair { //TODO rename into something like interval
  #ifdef SERIARIALIZATION_ENABLED
    friend class icecube::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
  #endif //SERIARIALIZATION_ENABLED  
  public: //properties
    ///down a string
    double minus_;
    ///up a string
    double plus_;
    
  private:
    ///blank constructor
    LimitPair();
      
  public:
    /// constructor
    /// @param m distance down a string
    /// @param p distance up a string
    LimitPair(const double m, const double p);
    
    ///the value is within the specified ring limits
    /// \param val the value/distance; a (neg) double value
    bool Within(const double val) const;
  };

  //=================== CLASS RingLimits ==============
  /// collection of limit-pairs as extending to further out rings
  class RingLimits {
  #ifdef SERIARIALIZATION_ENABLED
    friend class icecube::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
  #endif //SERIARIALIZATION_ENABLED
  public://TODO move private
    ///holds the limit pairs
    std::vector<LimitPair> limitPairs_;
  public:
    ///blank constructor
    RingLimits();
    
    ///construct directly
    RingLimits(const std::vector<LimitPair> &l);
    
    ///append a LimitPair
    void AddLimitPair(const LimitPair &lp); //FIXME away or make specific, so to add pair at position
    
    /** @brief How many Rings are configured 
      * @return -1 if there is nothing configured, 
      *          0 if only the central string is configured
      *          int>0 for number of configured Rings
      */
    int NRings() const;
    
    /// get the limit pair on a specific ring
    LimitPair GetLimitsOnRing(const int r) const;
  };
  
public:
  /// PARAM: holding the Hive Geometry
  const hive::HiveTopologyConstPtr hivetopo_;
  
  /// PARAM: specifies all available DOMs to connect from; object is function of signature: bool (const OMKey&)
  boost::function<bool (const OMKey&)> connectFrom_;
  /// PARAM: specifies all available DOMs to connect to; object is function of signature: bool (const OMKey&)
  boost::function<bool (const OMKey&)> connectTo_;
  /// PARAM: limits for the single, double and tripple dense light connection
  RingLimits ringLimits_;
  /// PARAM: Allow for DOMs to self-connect
  bool selfconnect_;
  /// PARAM: Connect these DOMs mutually; A->B => B->A 
  bool mutuallyconnect_;

public: //interface
  /// Constructor 
  HiveRelationConfig(const hive::HiveTopologyConstPtr hivetopo);
  
  /** @brief Build the DistanceMap and VicinityMap from the parameters
   * keep things simple and call this function from Geometry-frame
   * \param geo Detector geometry
   * \param calib Detector Calibration
   * \param status Detector Status
   * \param hashService the hasher for OMkeys; is created if pointer not set
   * \param distService the distance service; is creater if pointer not set
   */
  RelationPtr BuildRelation (
    const HashedGeometryConstPtr& hashedGeo) const;
};

typedef boost::shared_ptr<HiveRelationConfig> HiveRelationConfigPtr;
typedef boost::shared_ptr<const HiveRelationConfig> HiveRelationConfigConstPtr;


//===============================================================
//===================== IMPLEMENTATION ==========================
//===============================================================


//====================== CLASS LimitPair ========================
#ifdef SERIARIALIZATION_ENABLED
template <class Archive>
void HiveRelationConfig::LimitPair::serialize(Archive & ar, const unsigned int version) {
  ar & SERIALIZATION_MAKE_NVP(minus_);
  ar & SERIALIZATION_MAKE_NVP(plus_);
};
#endif //SERIARIALIZATION_ENABLED

inline
bool HiveRelationConfig::LimitPair::Within(const double val) const
  {return minus_ <= val && val <= plus_;};

  
//====================== CLASS RingLimits ========================
#ifdef SERIARIALIZATION_ENABLED
template<class Archive>
void HiveRelationConfig::RingLimits::serialize(Archive & ar, const unsigned int version) {
  ar & SERIALIZATION_MAKE_NVP(limitPairs_);
};
#endif //SERIARIALIZATION_ENABLED

#endif //RELATIONCONFIG_H
