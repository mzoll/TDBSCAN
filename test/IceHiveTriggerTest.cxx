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

#include "IceHiveZ/modules/IceHiveTrigger.h"

#include <boost/make_shared.hpp>
#include <boost/foreach.hpp>

#include "TestHelpers.h"
#include "ToolZ/IC86Topology.h"

using namespace std;
using namespace boost;
using namespace HitSorting;

using namespace hivetrigger;

TEST_GROUP(IceHiveTrigger);

const I3GeometryConstPtr geometry = boost::make_shared<const I3Geometry>(IC86Topology::Build_IC86_Geometry());
const I3CalibrationConstPtr calibration(new I3Calibration());
const I3DetectorStatusConstPtr status(new I3DetectorStatus());

const HashedGeometryConstPtr hashedGeo = boost::make_shared<const HashedGeometry>(geometry->omgeo);

TEST(IceHiveTrigger) {
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

  IceHiveTrigger iht(ht_param_set);
  
  I3DOMLaunchSeriesMap launchMap(GenerateDetectorNoiseDOMLaunches(10*I3Units::ms));
  
  typedef std::list<I3DOMLaunch_HitObject> I3DOMLaunch_HitObjectList;
  I3DOMLaunch_HitObjectList hitobj = OMKeyMap_To_HitObjects<I3DOMLaunch, I3DOMLaunch_HitObjectList>(launchMap);
  
  BOOST_FOREACH(const I3DOMLaunch_HitObject& ho, hitobj)
    iht.Feed(ho.GetOMKey(), ho.GetResponseObj(), ho.GetDAQTicks());

  iht.AdvanceUntil(std::numeric_limits<DAQTicks>::max());
  ENSURE_EQUAL(iht.FinalizedUntil(), std::numeric_limits<DAQTicks>::max());
  
  iht.GetTriggers();
};
