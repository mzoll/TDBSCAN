/**
 * $Id: HiveCleaningTest.cxx 153794 2017-03-09 10:43:19Z mzoll $
 * $Author: mzoll $
 * $Date: 2017-03-09 11:43:19 +0100 (Thu, 09 Mar 2017) $
 * $Revision: 153794 $
 *
 * A Unit test which generates some artificial test cases and let the Cleaning gnaw on them;
 */

#include <I3Test.h>

#include "IceHiveZ/algorithms/HiveCleaning.h"

#include "ToolZ/IC86Topology.h"

#include "TestHelpers.h"
#include "dataclasses/I3Constants.h"
#include "icetray/I3Units.h"

#include <boost/make_shared.hpp>
#include <boost/foreach.hpp>

using namespace std;
using namespace HitSorting;

TEST_GROUP(HiveCleaning);

const I3GeometryConstPtr geometry = boost::make_shared<I3Geometry>(IC86Topology::Build_IC86_Geometry());
const I3CalibrationConstPtr calibration(new I3Calibration());
const I3DetectorStatusConstPtr status(new I3DetectorStatus());

const HashedGeometryConstPtr hashedGeo = boost::make_shared<const HashedGeometry>(geometry->omgeo);

TEST(HiveCleaning) {
  
  /* ALTERNATIVE I
  //prepare a frame containing pulses
  I3FramePtr frame = boost::make_shared<I3Frame>(I3Frame::Physics);
  I3RecoPulseSeriesMapPtr recomap = boost::make_shared<I3RecoPulseSeriesMap>(GenerateTestRecoPulses());
  frame.Put("Pulses", recomap);

  //obtain hits
  OMKeyMap_HitFacility<I3RecoPulse> hitFacility(hasher, frame, "Pulses");
  const HitSet hits = hitFacility.GetHits<AbsHitSet>();
  */
  
  /* ALTERNATIVE II */
  I3RecoPulseSeriesMap recoMap = GenerateDetectorNoiseRecoPulses(1E6);
  
  typedef std::list<HitObject<I3RecoPulse> > I3RecoPulseHitObjectList;
  I3RecoPulseHitObjectList hol = OMKeyMap_To_HitObjects<I3RecoPulse, I3RecoPulseHitObjectList>(recoMap);
  
  const HitSet hits = HitObjects_To_Hits<I3RecoPulseHitObjectList, HitSet>(hol, hashedGeo->GetHashService());
  
  //prepare the settings for the HiveSplitter
  HiveCleaning_ParameterSet hc_param_set;
//   hc_param_set.multiplicity;
//   hc_param_set.max_tresidual_early;
//   hc_param_set.max_tresidual_late;
  hc_param_set.connectorBlock = boost::make_shared<ConnectorBlock>(hashedGeo);
  hc_param_set.connectorBlock->AddConnector(boost::make_shared<Connector>("ConnectAll",
                                                                          hashedGeo,
                                                                          boost::make_shared<BoolConnection>(hashedGeo, true),
                                                                          boost::make_shared<Relation>(hashedGeo->GetHashService(), true)));

  HiveCleaning hiveCleaning( hc_param_set );
  
  //perform the cleaning
  AbsHitSet cleanHits = hiveCleaning.Clean(hits);

  //everything should be disconnected, so no hits written out
  ENSURE_EQUAL(cleanHits.size(), hits.size(), "Cleaned Series has the same size, as nothing should be cleaned away");
};
