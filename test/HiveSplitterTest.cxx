/**
 * \file HiveSplitterTest.cxx
 * 
 * $Id: HiveSplitterTest.cxx 153794 2017-03-09 10:43:19Z mzoll $
 * $Author: mzoll $
 * $Date: 2017-03-09 11:43:19 +0100 (Thu, 09 Mar 2017) $
 * $Revision: 153794 $
 *
 * A Unit test which generates some artificial test cases and let the Splitter gnaw on them;
 */

#include <I3Test.h>

#include "IceHiveZ/algorithms/HiveSplitter.h"

#include "ToolZ/HitSorting.h"
#include "ToolZ/IC86Topology.h"

#include "dataclasses/I3Constants.h"
#include "icetray/I3Units.h"

#include <boost/make_shared.hpp>
#include <boost/foreach.hpp>

#include "TestHelpers.h"

using namespace std;
using namespace boost;
using namespace HitSorting;

TEST_GROUP(HiveSplitter);

//prepare some faky input objects
const I3GeometryConstPtr geometry = boost::make_shared<I3Geometry>(IC86Topology::Build_IC86_Geometry());
const I3CalibrationConstPtr calibration(new I3Calibration());
const I3DetectorStatusConstPtr status(new I3DetectorStatus());

const HashedGeometryConstPtr hashedGeo = boost::make_shared<const HashedGeometry>(geometry->omgeo);

TEST(HiveSplitting) {
  //perform one round splitting


  
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
  I3RecoPulseSeriesMap recoMap = GenerateDetectorNoiseRecoPulses(10*I3Units::ms);
  
  typedef std::list<HitObject<I3RecoPulse> > I3RecoPulseHitObjectList;
  I3RecoPulseHitObjectList hol = OMKeyMap_To_HitObjects<I3RecoPulse, I3RecoPulseHitObjectList>(recoMap);
  
  const HitSet hits = HitObjects_To_Hits<I3RecoPulseHitObjectList, HitSet>(hol, hashedGeo->GetHashService());

  //prepare the settings for the HiveSplitter
  hivesplitter::HiveSplitter_ParameterSet hs_param_set;
//   hs_param_set.multiplicity;
//   hs_param_set.multiplicityTimeWindow;
//   hs_param_set.acceptTimeWindow;
//   hs_param_set.rejectTimeWindow;
//   hs_param_set.mergeOverlap;
  hs_param_set.connectorBlock = boost::make_shared<ConnectorBlock>(hashedGeo);
  hs_param_set.connectorBlock->AddConnector(boost::make_shared<Connector>("ConnectNone",
                                                                          hashedGeo,
                                                                          boost::make_shared<BoolConnection>(hashedGeo, false),
                                                                          boost::make_shared<Relation>(hashedGeo->GetHashService(), false)));
  HiveSplitter hiveSplitter_nonecon( hs_param_set );

  //perform the splitting
  AbsHitSetSequence subEvents_nonecon = hiveSplitter_nonecon.Split(hits);

  //everything should be disconnected, so no hits written out
  ENSURE_EQUAL(subEvents_nonecon.size(), 0, "No SubEvents written out, as hits were not connected");
  
  
  //make a new round, this time everything is connected
  hs_param_set.connectorBlock = boost::make_shared<ConnectorBlock>(hashedGeo);
  hs_param_set.connectorBlock->AddConnector(boost::make_shared<Connector>("ConnectAll",
                                                                          hashedGeo,
                                                                          boost::make_shared<BoolConnection>(hashedGeo, true),
                                                                          boost::make_shared<Relation>(hashedGeo->GetHashService(), true)));  
  HiveSplitter hiveSplitter_allcon( hs_param_set );
  
  //perform the splitting
  AbsHitSetSequence subEvents_allcon = hiveSplitter_allcon.Split(hits);

  //everything should be disconnected, so no hits written out
//   ENSURE_EQUAL(subEvents_allcon.size(), 1, "One SubEvents written out, as hits were all connected");
//   ENSURE_EQUAL(subEvents_allcon.begin()->size(), hits.size(), "All Hits are contained in this single subevent");
};
