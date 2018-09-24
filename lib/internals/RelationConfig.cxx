/**
 * \file HiveRelationConfig.cxx
 *
 * (c) 2012 the IceCube Collaboration
 *
 * $Id: HiveSplitter.h 99900 2013-02-26 10:10:43Z mzoll $
 * \version $Revision: 99900 $
 * \date $Date: 2013-02-26 11:10:43 +0100 (Tue, 26 Feb 2013) $
 * \author Marcel Zoll <marcel.zoll@fysik.su.se>
 *
 */

#include "IceHiveZ/internals/RelationConfig.h"

#include <cmath>

#define MULTITHREADING_ENABLED 0

#if MULTITHREADING_ENABLED
  #include <thread>         // std::thread
  #include <mutex>          // std::mutex

  static std::mutex mtx;
#endif //MULTITHREADING_ENABLED

using namespace std;

//===================== CLASS RelationConfig ==================

//================ CLASS SimpleRelationConfig ==================

SimpleRelationConfig::SimpleRelationConfig(
  boost::function<bool (const OMKey&, const OMKey&)> callobj)
: callobj_(callobj)
{};

RelationPtr SimpleRelationConfig::BuildRelation (
  const HashedGeometryConstPtr& hashedGeo) const
{ return boost::make_shared<Relation>(hashedGeo->GetHashService(), callobj_); };


//===================== STRUCT LimitPairs ==================

HiveRelationConfig::LimitPair::LimitPair() {};
    
HiveRelationConfig::LimitPair::LimitPair(const double m, const double p):
  minus_(m), plus_(p)
{};


//=================== CLASS RingLimits ==============

HiveRelationConfig::RingLimits::RingLimits():
  limitPairs_()
{};

///construct directly
HiveRelationConfig::RingLimits::RingLimits(const std::vector<LimitPair> &l):
  limitPairs_(l)
{};

void HiveRelationConfig::RingLimits::AddLimitPair(const LimitPair &lp) {
  limitPairs_.push_back(lp);
};

int HiveRelationConfig::RingLimits::NRings() const
  {return limitPairs_.size()-1;};
    
HiveRelationConfig::LimitPair HiveRelationConfig::RingLimits::GetLimitsOnRing(const int r) const 
  {return limitPairs_.at(r);};
  

//=============CLASS HiveRelationConfig ========================

HiveRelationConfig::HiveRelationConfig(
  const hive::HiveTopologyConstPtr hivetopo)
: hivetopo_(hivetopo)
{};
  
RelationPtr HiveRelationConfig::BuildRelation (
  const HashedGeometryConstPtr& hashedGeo) const   
{
  //verify helper objects
  if (!hivetopo_)
    log_fatal("No Hive set");
    
  const CompactOMKeyHashServiceConstPtr hasher = hashedGeo->GetHashService();
  const PositionServiceConstPtr posService = hashedGeo->GetPosService();
  
  //Start the map creation process
  RelationPtr rs = boost::make_shared<Relation>(hasher, false);
  
  log_info("Constructing Relation from HiveRelationConfig");
  //Construct the relation map:
  // For this fill the boolmap with all connections that can be made, but never reset them!
  
  using namespace hive;
  
  //NOTE as the following loops are static, they can in principle be multithreaded
  // however the access and manipulation of the result-type (Relation->IndexMap->boost::dynamic_bitset)
  // might not neccessarily be thread-save
  
  //===== LOOP A: ConnectFrom =====
  for (unsigned matrix_x=0; matrix_x<hasher->HashSize(); ++matrix_x) {
    const OMKey omkey_A = hasher->OMKeyFromHash(matrix_x);
    log_debug_stream("====Fill next row === : "<<matrix_x<<":"<<omkey_A);

    if (! connectFrom_(omkey_A))
      continue;

    const unsigned center_string = omkey_A.GetString();
    if (! hivetopo_->HoldsCenterString(center_string))
      continue;

    const I3Position& pos_A = posService->GetPosition(matrix_x);
    const double z_A(pos_A.GetZ());

    //===== LOOP B : ConnectTo=====
    for (unsigned matrix_y=0; matrix_y<hasher->HashSize(); ++matrix_y) {
      const OMKey omkey_B = hasher->OMKeyFromHash(matrix_y);
      log_debug_stream("Looking for " << omkey_A << " and " << omkey_B);
      
      if (selfconnect_ && matrix_x==matrix_y) {
        rs->SetRelated(matrix_x, matrix_y, true);
        continue;
      }
      
      if (! connectTo_(omkey_B))
        continue;
       
      const unsigned lookup_string = omkey_B.GetString();

      //Analyse the Ring relation between the srings of the two DOMs
      const int ring = hivetopo_->WhichRing(center_string, lookup_string);
      if (ring == -1) {//not in the ring indexing range
        log_trace("Not included in ring index range");
        continue;
      }
      log_trace_stream("Looking for " << omkey_A << " and " << omkey_B << " ring "<<ring << " max"<<ringLimits_.NRings());

      if (ring > ringLimits_.NRings()) {
        log_trace_stream("Ring "<<ring<<" too far away; max Ringlimits "<<ringLimits_.NRings());
        continue;
      }
      
      const I3Position& pos_B = posService->GetPosition(matrix_y);
      const double z_B(pos_B.GetZ());
        
      const double& minusLimit = ringLimits_.GetLimitsOnRing(ring).minus_;
      const double& plusLimit = ringLimits_.GetLimitsOnRing(ring).plus_;
      log_trace_stream("ring "<< ring << " z_A "<< z_A << " z_B " << z_B << " ringlimit- " <<minusLimit<< " ringlimit+ " <<plusLimit);

      const double zdist_AB = z_B-z_A;
      if (std::isnan(minusLimit) && std::isnan(plusLimit)) {
        log_trace("Not configured rings");
      }
      else if (minusLimit<= zdist_AB && zdist_AB <= plusLimit) {
        rs->SetRelated(matrix_x, matrix_y, true);
        if (mutuallyconnect_)
          rs->SetRelated(matrix_y, matrix_x, true);
        log_trace("DOMs are connected");
      }
      else {
        log_trace("Not included in ring limits");
      }
    } //LOOP_B
  } //LOOP_A
  log_info("DONE Constructing RelationMap");
  log_debug("Leaving BuildDistanceMap()");
  
  return rs;
};
