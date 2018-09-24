/**
 * $Id: TestHelpers.cxx 124115 2014-10-03 16:42:30Z mzoll $
 * $Author: mzoll $
 * $Date: 2014-10-03 18:42:30 +0200 (fre, 03 okt 2014) $
 * $Revision: 124115 $
 *
 * A Unit test which generates some artificial test cases and let the Cleaning gnaw on them;
 */

#define ICINTERSTRINGSPACING 125.
#define DCINTERSTRINGSPACING 72.17

#include "dataclasses/I3Constants.h"
#include "icetray/I3Units.h"

#include "dataclasses/physics/I3RecoPulse.h"
#include "dataclasses/physics/I3DOMLaunch.h"

#include "IceHiveZ/__SERIALIZATION.h"

///directly init all relevant fields of a I3RecoPulse
I3RecoPulse MakeRecoPulse (const double t, const double c, const double w=1., const uint8_t flags=0);

///generate a number of RecoPulses in a I3RecoPulseSeriesMap
I3RecoPulseSeriesMap GenerateTestRecoPulses();

///generate a number of RecoPulses in a I3RecoPulseSeriesMap of detector Noise
I3RecoPulseSeriesMap GenerateDetectorNoiseRecoPulses(const double max_time_range);

///directly init all relevant fields of a I3DOMLaunch
I3DOMLaunch MakeDOMLaunch (const double t, const bool lcbit=false);

///generate a number of DOMLaunches in a I3DOMLaunchSeriesMap
I3DOMLaunchSeriesMap GenerateTestDOMLaunches();

///generate a number of RecoPulses in a I3RecoPulseSeriesMap of detector Noise
I3DOMLaunchSeriesMap GenerateDetectorNoiseDOMLaunches(const double max_time_range);


#if SERIALIZATION_ENABLED
template <class T>
void serialize_object(
  const T* t_save,
  T*& t_load)
{
  std::stringstream ss;
  SERIALIZATION_NS_BASE::archive::portable_binary_oarchive oa(ss);
//   SERIALIZATION_NS_BASE::archive::xml_oarchive oa(ss);
  oa << t_save;
  SERIALIZATION_NS_BASE::archive::portable_binary_iarchive ia(ss);
//   SERIALIZATION_NS_BASE::archive::xml_iarchive ia(ss);
  ia >> t_load;
};

template <class T>
void serialize_object(
  const boost::shared_ptr<T> t_save,
  boost::shared_ptr<T>& t_load)
{
  std::stringstream ss;
  SERIALIZATION_NS_BASE::archive::portable_binary_oarchive oa(ss);
//   SERIALIZATION_NS_BASE::archive::xml_oarchive oa(ss);
  oa << t_save;
  SERIALIZATION_NS_BASE::archive::portable_binary_iarchive ia(ss);
//   SERIALIZATION_NS_BASE::archive::xml_iarchive ia(ss);
  ia >> t_load;
};
#endif //SERIALIZATION_ENABLED


#include "ToolZ/OMKeyHash.h"
CompactOMKeyHashServiceConstPtr
DummyHashService(const size_t size);


