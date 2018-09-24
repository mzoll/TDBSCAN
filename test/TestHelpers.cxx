/**
 * $Id: TestHelpers.cxx 124115 2014-10-03 16:42:30Z mzoll $
 * $Author: mzoll $
 * $Date: 2014-10-03 18:42:30 +0200 (fre, 03 okt 2014) $
 * $Revision: 124115 $
 *
 * A Unit test which generates some artificial test cases and let the Cleaning gnaw on them;
 */

#include "TestHelpers.h"

#include <gsl/gsl_rng.h>

#include "boost/make_shared.hpp"

#include "ToolZ/OMKeyHash.h"
#include "ToolZ/HitSorting.h"


I3RecoPulse MakeRecoPulse (const double t, const double c, const double w, const uint8_t flags) {
  I3RecoPulse p;
  p.SetTime(t);
  p.SetCharge(c);
  p.SetWidth(w);
  p.SetFlags(flags);
  return p;
}

I3RecoPulseSeriesMap GenerateTestRecoPulses(
  const size_t n_gen,
  const double timerange
) {
  I3RecoPulseSeriesMap pulseMap;
  for (uint64_t i=1; i<=86; i++) {
    for (uint64_t j=1; j<=60; j++) {
      for (uint64_t k=0; k<20; k++) {
        pulseMap[OMKey(i,j)].push_back(MakeRecoPulse(k, k));
      }
    }
  }
  return pulseMap;
};

I3RecoPulseSeriesMap GenerateDetectorNoiseRecoPulses(const double max_time_range_ns) {
  using namespace HitSorting;
  
  assert(max_time_range_ns>0.);
  
  const double DOM_NOISE_RATE_IC = 500.; //Hz
  const double DOM_NOISE_RATE_DC = 800.; //Hz
  
  const uint64_t N_HITS_EXPECTED_IC = uint64_t(max_time_range_ns*1.E-9*DOM_NOISE_RATE_IC*78*60);
  const uint64_t N_HITS_EXPECTED_DC = uint64_t(max_time_range_ns*1.E-9*DOM_NOISE_RATE_DC*8*60);
  
  //create hashable range
  std::set<OMKey> omkey_set;
  for (uint64_t i=1; i<=86; i++) {
    for (uint64_t j=1; j<=60; j++) {
      omkey_set.insert(omkey_set.end(), OMKey(i,j));
    }
  }
  const CompactOMKeyHashServiceConstPtr hasher
    = boost::make_shared<CompactOMKeyHashService>(omkey_set);
  
  //initialize the random number generator
  const gsl_rng_type * T;
  gsl_rng * r;
  gsl_rng_env_setup();
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);

  //typedef std::list<HitObject<I3RecoPulse> > HitObjectsList;
  typedef std::list<HitObjectOriginal<I3RecoPulse> > HitObjectsList;
  HitSet hits;
  
  //for IceCube DOMs
  HitObjectsList hitobjects_ic;
  while (hitobjects_ic.size()<N_HITS_EXPECTED_IC) {
    uint64_t string_ic = gsl_rng_uniform_int (r, 77)+1;
    uint64_t om = gsl_rng_uniform_int (r, 59)+1;
    uint64_t time = gsl_rng_uniform (r) * max_time_range_ns;
    
    hitobjects_ic.push_back(HitObjectOriginal<I3RecoPulse>(OMKey(string_ic,om),MakeRecoPulse(time, 1)));
    hits.insert(Hit(hasher->HashFromOMKey(hitobjects_ic.rbegin()->GetOMKey()), time, hitobjects_ic.back()));
  }
  
  //for DeepCore DOMs
  HitObjectsList hitobjects_dc;
  while (hitobjects_dc.size()<N_HITS_EXPECTED_DC) {
    uint64_t string_dc = gsl_rng_uniform_int (r, 8)+79;
    uint64_t om = gsl_rng_uniform_int (r, 59)+1;
    uint64_t time = gsl_rng_uniform (r) * max_time_range_ns;
    
    hitobjects_dc.push_back(HitObjectOriginal<I3RecoPulse>(OMKey(string_dc,om),MakeRecoPulse(time, 1)));    
    hits.insert(Hit(hasher->HashFromOMKey(hitobjects_dc.rbegin()->GetOMKey()), time, hitobjects_dc.back()));
  }
  
  gsl_rng_free (r);
  log_info_stream("Delivered "<<hits.size()<<" detector noise hits");
  
  I3RecoPulseSeriesMap recoMap = HitSorting::Hits_To_OMKeyMap<I3RecoPulse, HitSet>(hits, hasher);

  return recoMap;
};

I3DOMLaunch MakeDOMLaunch (const double t, const bool lcbit)
{
  I3DOMLaunch l;
  l.SetStartTime(t);
  l.SetLCBit(lcbit);
//     void SetStartTime(double starttime)
//     void SetTriggerType(TriggerType trigger)
//     void SetTriggerMode(TriggerMode situation)
//     void SetWhichATWD(ATWDselect WhichATWD)
//     void SetRawATWD( std::vector<std::vector<int> >& v )
//     void SetRawFADC( std::vector<int>& v)
//     void SetLCBit(bool LCBit)
//     void SetIsPedestalSub(bool Pedestal)
//     void SetChargeStampHighestSample(unsigned int highsample)
//     void SetRawATWDChargeStamp(unsigned int chargesum)
//     void SetWhichATWDChargeStamp(unsigned int channelid)
  return l;
}

I3DOMLaunchSeriesMap GenerateTestDOMLaunches() {
  I3DOMLaunchSeriesMap launchMap;
  for (uint64_t i=1; i<=86; i++) {
    for (uint64_t j=1; j<=60; j++) {
      for (uint64_t k=0; k<20; k++) {
        launchMap[OMKey(i,j)].push_back(MakeDOMLaunch(k));
      }
    }
  }
  return launchMap;
};

I3DOMLaunchSeriesMap GenerateDetectorNoiseDOMLaunches(const double max_time_range_ns) {
  using namespace HitSorting;
  
  assert(max_time_range_ns>0.);
  
  const double DOM_NOISE_RATE_IC = 500.; //Hz
  const double DOM_NOISE_RATE_DC = 800.; //Hz
  
  const uint64_t N_HITS_EXPECTED_IC = uint64_t(max_time_range_ns*1.E-9*DOM_NOISE_RATE_IC*78*60);
  const uint64_t N_HITS_EXPECTED_DC = uint64_t(max_time_range_ns*1.E-9*DOM_NOISE_RATE_DC*8*60);
  
  //create hashable range
  std::set<OMKey> omkey_set;
  for (uint64_t i=1; i<=86; i++) {
    for (uint64_t j=1; j<=60; j++) {
      omkey_set.insert(omkey_set.end(), OMKey(i,j));
    }
  }
  const CompactOMKeyHashServiceConstPtr hasher
    = boost::make_shared<CompactOMKeyHashService>(omkey_set);
  
  //initialize the random number generator
  const gsl_rng_type * T;
  gsl_rng * r;
  gsl_rng_env_setup();
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);

  //typedef std::list<HitObject<I3DOMLaunch> > HitObjectsList;
  typedef std::list<HitObjectOriginal<I3DOMLaunch> > HitObjectsList;
  HitSet hits;
  
  //for IceCube DOMs
  HitObjectsList hitobjects_ic;
  while (hitobjects_ic.size()<N_HITS_EXPECTED_IC) {
    uint64_t string_ic = gsl_rng_uniform_int (r, 77)+1;
    uint64_t om = gsl_rng_uniform_int (r, 59)+1;
    uint64_t time = gsl_rng_uniform (r) * max_time_range_ns;
    
    hitobjects_ic.push_back(HitObjectOriginal<I3DOMLaunch>(OMKey(string_ic,om),MakeDOMLaunch(time)));
    hits.insert(Hit(hasher->HashFromOMKey(hitobjects_ic.rbegin()->GetOMKey()), time, hitobjects_ic.back()));
  }
  
  //for DeepCore DOMs
  HitObjectsList hitobjects_dc;
  while (hitobjects_dc.size()<N_HITS_EXPECTED_DC) {
    uint64_t string_dc = gsl_rng_uniform_int (r, 8)+79;
    uint64_t om = gsl_rng_uniform_int (r, 59)+1;
    uint64_t time = gsl_rng_uniform (r) * max_time_range_ns;
    
    hitobjects_dc.push_back(HitObjectOriginal<I3DOMLaunch>(OMKey(string_dc,om),MakeDOMLaunch(time)));    
    hits.insert(Hit(hasher->HashFromOMKey(hitobjects_dc.rbegin()->GetOMKey()), time, hitobjects_dc.back()));
  }
  
  gsl_rng_free (r);
  log_info_stream("Delivered "<<hits.size()<<" detector noise hits");
  
  I3DOMLaunchSeriesMap recoMap = HitSorting::Hits_To_OMKeyMap<I3DOMLaunch, HitSet>(hits, hasher);

  return recoMap;
};


//create some (global) objects which are used in all the hashings
// I3GeometryConstPtr geo = boost::make_shared<const I3Geometry>(IC86Topology::Build_IC86_Geometry());
// CompactOMKeyHashServiceConstPtr hasher = boost::make_shared<CompactOMKeyHashService>(geo);
// DistanceServiceConstPtr distService = boost::make_shared<DistanceService>(geo,hasher); 

CompactOMKeyHashServiceConstPtr
DummyHashService(const size_t size) {
  std::set<OMKey> omkey_set;
  for (size_t i=0; i<size; i++)
    omkey_set.insert(omkey_set.end(), OMKey(i,i));
  
  return boost::make_shared<const CompactOMKeyHashService>(omkey_set);
};
