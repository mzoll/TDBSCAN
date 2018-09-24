/**
 * \file RelationTest.cxx
 *
 * (c) 2013 the IceCube Collaboration
 *
 * $Id: OMKeyHashTest.cxx 148312 2016-07-11 17:41:28Z mzoll $
 * \version $Revision: 148312 $
 * \date $Date: 2016-07-11 19:41:28 +0200 (mån, 11 jul 2016) $
 * \author mzoll <marcel.zoll@fysik.su.se>
 *
 * Unit test to test the 
 */

#include <I3Test.h>

#include "IceHiveZ/internals/Relation.h"

#include "ToolZ/IC86Topology.h"

#include "TestHelpers.h"

TEST_GROUP(Relation);

//prepare some faky input objects
// static I3GeometryConstPtr geometry = boost::make_shared<const I3Geometry>(IC86Topology::Build_IC86_Geometry());
// static I3CalibrationConstPtr calibration(new I3Calibration());
// static I3DetectorStatusConstPtr status(new I3DetectorStatus());
// static CompactOMKeyHashServiceConstPtr hashService = boost::make_shared<const CompactOMKeyHashService>(geometry);
static CompactOMKeyHashServiceConstPtr hashService =  DummyHashService(100);

TEST(Construct_simple) {
  Relation rel_none(hashService, false);
  
  //just test the first and last entry
  const OMKey minOMKey = hashService->OMKeyFromHash(0);
  const OMKey maxOMKey = hashService->OMKeyFromHash(hashService->HashSize()-1);
  ENSURE_EQUAL(rel_none.AreRelated(minOMKey, minOMKey), false);
  ENSURE_EQUAL(rel_none.AreRelated(maxOMKey, maxOMKey), false);
  
  Relation rel_all(hashService, true);
  
  //just test the first and last entry
  ENSURE_EQUAL(rel_all.AreRelated(minOMKey, minOMKey), true);
  ENSURE_EQUAL(rel_all.AreRelated(maxOMKey, maxOMKey), true);
};


TEST (Construct_Predicate) {
  //make a struct that evaluates true if both OMKeys are even OM numbers
  struct SetEven {
    bool operator() (const OMKey& omkey_A, const OMKey& omkey_B) const 
      {return (!(omkey_A.GetOM()%2) && !(omkey_B.GetOM()%2));};
  };
  SetEven se;
  
  Relation rel(hashService, se);
  
  for (uint64_t i=0; i<hashService->HashSize(); i++) {
    const OMKey omkey_A= hashService->OMKeyFromHash(i);
    for (uint64_t j=0; j<hashService->HashSize(); j++) {
      const OMKey omkey_B = hashService->OMKeyFromHash(j);
      ENSURE_EQUAL(rel.AreRelated(omkey_A, omkey_B), se(omkey_A, omkey_B));
    }
  }
};


TEST(SetAndGet) {
  Relation rel(hashService, false);
  const OMKey minOMKey = hashService->OMKeyFromHash(0);
  const OMKey maxOMKey = hashService->OMKeyFromHash(hashService->HashSize()-1);
  
  rel.SetRelated(minOMKey, minOMKey, true);
  rel.SetRelated(maxOMKey, maxOMKey, true);
  
  ENSURE_EQUAL(rel.AreRelated(minOMKey, minOMKey), true);
  ENSURE_EQUAL(rel.AreRelated(maxOMKey, maxOMKey), true);
}


#if SERIALIZATION_ENABLED
TEST(Serialize_raw_ptr){
  Relation* rel_save = new Relation(hashService, false);
  Relation* rel_load = NULL;
  
  serialize_object(rel_save, rel_load);
  delete rel_load;
};

TEST(Serialize_boost_shared_ptr){
  RelationPtr rel_save= boost::make_shared<Relation>(hashService, false);
  RelationPtr rel_load;
  
  serialize_object(rel_save, rel_load);
};
#endif //SERIALIZATION_ENABLED
