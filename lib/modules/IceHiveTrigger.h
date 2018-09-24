/**
 * \file IceHiveTrigger.h
 *
 * (c) 2012 the IceCube Collaboration
 *
 * $Id: IceHiveTrigger.h 99900 2013-02-26 10:10:43Z mzoll $
 * \version $Revision: 99900 $
 * \date $Date: 2013-02-26 11:10:43 +0100 (Tue, 26 Feb 2013) $
 * \author Marcel Zoll <marcel.zoll@fysik.su.se>
 *
 * A interface between the HiveTrigger algorithm and the LaunchStreamTrigger interface
 */

#ifndef ICEHIVETRIGGER_H
#define ICEHIVETRIGGER_H

#include <boost/make_shared.hpp>

#include "IceHiveZ/algorithms/HiveTrigger.h"

#include "ToolZ/I3RUsageTimer.h"

#include "hitspool-reader/HitSpoolTrigger.h"

///The main module which unites the algorithms HiveSplitter and TriggerSplitter
class IceHiveTrigger : public hitspooltrigger::LaunchStreamTrigger {
  SET_LOGGER("IceHiveTrigger");

protected://parameters
  //========================
  // Configurable Parameters
  //========================  
  ///PARAM: which are delivered to HiveSplitter
  hivetrigger::HiveTrigger_ParameterSet ht_params_;
  /// PARAM: Minimal size on an subEvent to be considered as a Trigger
  size_t minEventSize_;
                  
private: //bookkeeping  
  //hits processed
  uint64_t n_hits_in_;
  //triggers produced
  uint64_t n_triggers_;
  //stopwatch
  I3RUsageTimer totRUsageEatTimer_;
  //stopwatch for splitter
  I3RUsageTimer totRUsageHiveTriggerTimer_;

private: //properties and methods related to configuration
  /// a global hasher to translate OMKeys to Hashes and vice versa
  CompactOMKeyHashServiceConstPtr hashService_;
  //facilitate the splitting
  ///most private HiveSplitter instance
  HiveTrigger* hiveTrigger_;
  /// holds all waiting triggerWindows; are merged and eventually pushed out
  hitspooltrigger::TriggerSet waitingTriggers_;
  
public: //methods
  //================
  // Main Interface
  //================
  /// Constructor: configure Default values, register Parameters, register Outbox
  IceHiveTrigger(
    const hivetrigger::HiveTrigger_ParameterSet& ht_params,
    const size_t minEventSize =1);
  /// Destructor
  virtual ~IceHiveTrigger();
  ///tell until which time the set of triggers is final
  hitspooltime::DAQTicks FinalizedUntil() const;
private:
  ///consumate the hit and do something with it
  void Eat(const OMKey& dom, const I3DOMLaunch& launch, const hitspooltime::DAQTicks DAQTime);
  ///collect all the triggers the deep deep hidden places
  void CollectTriggers();
  /// advance this to this time
  void AdvanceTime(const hitspooltime::DAQTicks DAQTime);
  /// report internal state of IceHive
  void ReportState() const;
};

typedef boost::shared_ptr<IceHiveTrigger> IceHiveTriggerPtr;
typedef boost::shared_ptr<const IceHiveTrigger> IceHiveTriggerConstPtr;

#endif
