/**
 * \file Hive-libTest.cxx
 *
 * (c) 2013 the IceCube Collaboration
 *
 * $Id: HiveTest.cxx 152659 2017-01-13 15:29:26Z mzoll $
 * \version $Revision: 152659 $
 * \date $Date: 2017-01-13 16:29:26 +0100 (Fri, 13 Jan 2017) $
 * \author mzoll <marcel.zoll@icecube.wisc.edu>
 *
 * This Unit test is to check some interactions of the Hive library. but nothing fancy
 */

#include <I3Test.h>

#include "IceHiveZ/internals/Hive.h"

#include "TestHelpers.h"

#include <boost/make_shared.hpp>

using namespace hive;

TEST_GROUP(Hive);

TEST(Build_a_hive) {
  
  StringRings sr(1); //create StringRings around string 36

  ENSURE_EQUAL(sr.GetCenter(), 1);
  ENSURE_EQUAL(sr.GetNRings(), 0);
  
  //create the first ring
  uint first_ring[] = {2};
  Ring ring1(first_ring, first_ring+1);
  sr.SetRing(1, ring1);
  ENSURE_EQUAL(sr.GetNRings(), 1);
  ENSURE_EQUAL(*sr.GetRing(0).begin(), 1);
  ENSURE_EQUAL(*sr.GetRing(1).begin(), 2);

  sr.AddStringToRing(3,2); //create second ring and add string 3 to it
  
  ENSURE_EQUAL(sr.GetNRings(), 2);
  ENSURE_EQUAL(*sr.GetRing(2).begin(), 3);

  ENSURE_EQUAL(sr.WhichRing(1), 0);
  ENSURE(sr.IsRingX(1,0));
  
  ENSURE_EQUAL(sr.WhichRing(2), 1);
  ENSURE(sr.IsRingX(2,1));
  
  ENSURE_EQUAL(sr.WhichRing(3), 2);
  ENSURE(sr.IsRingX(3,2));
  
  
  HiveTopology ht;
  
  ht.AddStringRing(sr);
  
  ENSURE(ht.HoldsCenterString(1));
  
  ENSURE_EQUAL(ht.WhichRing(1,1), 0);
  ENSURE(ht.IsRingX(1,1,0));
  
  ENSURE_EQUAL(ht.WhichRing(1,2), 1);
  ENSURE(ht.IsRingX(1,2,1));
  
  ENSURE_EQUAL(ht.WhichRing(1,3), 2);
  ENSURE(ht.IsRingX(1,3,2));

  StringRings sr4(4);
  
  ht.AddStringRing(sr4);
  ht.MutualAddStringToRing(1, 4, 3);
  
  ENSURE_EQUAL(ht.WhichRing(1,4), 3);
  ENSURE(ht.IsRingX(1,4,3));
  
  ENSURE_EQUAL(ht.WhichRing(4,1), 3);
  ENSURE(ht.IsRingX(4,1,3));
}

#if SERIALIZATION_ENABLED
TEST(StringRing_Serialize_raw_ptr){
  StringRings* ht_save = new StringRings(1);
  StringRings* ht_load = (StringRings*)malloc(sizeof(StringRings));
  
  serialize_object(ht_save, ht_load);
};

// TEST(StringRing_Serialize_boost_shared_ptr){
//   StringRingsPtr ht_save = boost::make_shared<StringRings>(1);
//   StringRingsPtr ht_load;
//   
//   serialize_object(ht_save, ht_load);
// };

TEST(HiveTopology_Serialize_raw_ptr){
  HiveTopology* ht_save = new HiveTopology();
  HiveTopology* ht_load = NULL;
  
  serialize_object(ht_save, ht_load);
  delete ht_save;
  delete ht_load;
};

TEST(HiveTopology_Serialize_boost_shared_ptr){
  HiveTopologyPtr ht_save = boost::make_shared<HiveTopology>();
  HiveTopologyPtr ht_load;
  
  serialize_object(ht_save, ht_load);
};
#endif //SERIALIZATION_ENABLED

