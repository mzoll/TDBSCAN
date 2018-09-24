/**
 * \file HiveTriggerTest.cxx
 * 
 * $Id: HiveTriggerTest.cxx 150050 2016-09-15 12:45:28Z mzoll $
 * $Author: mzoll $
 * $Date: 2016-09-15 14:45:28 +0200 (tor, 15 sep 2016) $
 * $Revision: 150050 $
 *
 * A Unit test which generates some artificial test cases and let the Splitter gnaw on them;
 */

#include <I3Test.h>

#include "IceHiveZ/algorithms/HiveTrigger.h"

#include <boost/make_shared.hpp>
#include <boost/foreach.hpp>

#include "TestHelpers.h"

#include "ToolZ/IC86Topology.h"

using namespace std;
using namespace boost;
using namespace HitSorting;

using namespace hivetrigger;

TEST_GROUP(HiveTrigger);

const I3GeometryConstPtr geometry = boost::make_shared<const I3Geometry>(IC86Topology::Build_IC86_Geometry());
const I3CalibrationConstPtr calibration(new I3Calibration());
const I3DetectorStatusConstPtr status(new I3DetectorStatus());

const HashedGeometryConstPtr hashedGeo = boost::make_shared<const HashedGeometry>(geometry->omgeo);

TEST(HiveTrigger) {
  //perform a round of triggering

  //prepare some faky input objects


  //prepare the settings for the HiveTrigger
  HiveTrigger_ParameterSet ht_param_set;
//   ht_param_set.multiplicity;
//   ht_param_set.multiplicityTimeWindow;
//   ht_param_set.acceptTimeWindow;
//   ht_param_set.rejectTimeWindow;
  ht_param_set.connectorBlock = boost::make_shared<ConnectorBlock>(hashedGeo);
  ht_param_set.connectorBlock->AddConnector(boost::make_shared<Connector>("ConnectNone",
                                                                          hashedGeo,
                                                                          boost::make_shared<BoolConnection>(hashedGeo, false),
                                                                          boost::make_shared<Relation>(hashedGeo->GetHashService(), true)));

  HiveTrigger hiveTrigger( ht_param_set );

  
  /* ALTERNATIVE I
  //prepare a frame containing pulses
  I3FramePtr frame = boost::make_shared<I3Frame>(I3Frame::Physics);
  I3DOMLaunchSeriesMapPtr launchmap = boost::make_shared<I3DOMLaunchSeriesMap>(GenerateTestDOMLaunches());
  frame.Put("Launches", recomap);

  //obtain hits
  OMKeyMap_HitFacility<I3DOMLaunch> hitFacility(hasher, frame, "Launches");
  const HitSorting::AbsDAQHitSet hits = hitFacility.GetAbsDAQHits<HitSorting::AbsDAQHitSet>();
  */
  
  /* ALTERNATIVE II */
  I3DOMLaunchSeriesMap launchMap(GenerateDetectorNoiseDOMLaunches(10*I3Units::ms));
  
  typedef std::list<I3DOMLaunch_HitObject> I3DOMLaunch_HitObjectList;
  const I3DOMLaunch_HitObjectList hitobjs = OMKeyMap_To_HitObjects<I3DOMLaunch, I3DOMLaunch_HitObjectList>(launchMap);
  
  const AbsDAQHitSet hits = HitObjects_To_AbsDAQHits<I3DOMLaunch_HitObjectList, AbsDAQHitSet>(hitobjs, hashedGeo->GetHashService());

  BOOST_FOREACH(const AbsDAQHit& h, hits)
    hiveTrigger.AddHit(h);

  hiveTrigger.FinalizeSubEvents();
  hiveTrigger.PullSubEvents();
};
