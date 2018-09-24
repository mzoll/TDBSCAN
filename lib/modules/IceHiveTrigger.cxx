/**
 * \file IceHiveTrigger.cxx
 *
 * (c) 2012 the IceCube Collaboration
 *
 * $Id: IceHiveTrigger.h 99900 2013-02-26 10:10:43Z mzoll $
 * \version $Revision: 99900 $
 * \date $Date: 2013-02-26 11:10:43 +0100 (Tue, 26 Feb 2013) $
 * \author Marcel Zoll <marcel.zoll@fysik.su.se>
 *
 * The IceTray I3Module wrapper around the central algorithm HiveSplitter and TriggerSplitter,
 * which in turn use a API trough SubEventStartStop
 */

#include "IceHiveZ/modules/IceHiveTrigger.h"

#include "icetray/I3Units.h"
#include "icetray/I3Int.h"
#include "dataclasses/I3MapOMKeyMask.h"
#include "dataclasses/physics/I3EventHeader.h"

using namespace hivetrigger;
using namespace HitSorting;

//===============class IceHiveTrigger=================================

IceHiveTrigger::IceHiveTrigger(
  const HiveTrigger_ParameterSet& ht_params,
  const size_t minEventSize):
  //parameter sets for subordinated modules
  ht_params_(ht_params),
  minEventSize_(0), 
  //bookkeeping
  n_hits_in_(0),
  n_triggers_(0),
  //initialize services
  hashService_(),
  //initialize the splitter algorithms
  hiveTrigger_(nullptr)
{
  log_debug("Creating IceHiveTrigger instance");
    //configuration needs to be done during init  
  hashService_ = ht_params_.connectorBlock->GetHashService();
  hiveTrigger_ = new HiveTrigger( ht_params_ );
  log_info("This is IceHiveTrigger!");
  
  log_debug("Leaving Init()");
}


IceHiveTrigger::~IceHiveTrigger() {
  if (hiveTrigger_!=NULL)
    delete hiveTrigger_;
  
  log_debug("Entering Finish()");
  
  log_notice_stream("Processed "<<n_hits_in_<<" hits producing "<<n_triggers_<<" triggers");
  
  I3RUsagePtr totalRUsage = totRUsageEatTimer_.GetTotalRUsage();
  log_notice(
      "%lu calls to Feed: %s, %.2fus per launch, %.2fms per trigger",
      n_hits_in_,
      convertI3RUsageToString(*totalRUsage).c_str(),
      (n_hits_in_ ? totalRUsage->wallclocktime/I3Units::microsecond/n_hits_in_ : 0),
      (n_triggers_ ? totalRUsage->wallclocktime/I3Units::millisecond/n_triggers_ : 0));
  
  totalRUsage = totRUsageHiveTriggerTimer_.GetTotalRUsage();
  log_notice(
      "%lu calls to HiveTrigger: %s, %.2fus per launch, %.2fms per trigger",
      n_hits_in_,
      convertI3RUsageToString(*totalRUsage).c_str(),
      (n_hits_in_ ? totalRUsage->wallclocktime/I3Units::microsecond/n_hits_in_ : 0),
      (n_triggers_ ? totalRUsage->wallclocktime/I3Units::millisecond/n_triggers_ : 0));
}

hitspooltime::DAQTicks IceHiveTrigger::FinalizedUntil() const
  {return hiveTrigger_->FinalizedUntil();};

void IceHiveTrigger::Eat(const OMKey& dom, const I3DOMLaunch& launch, const hitspooltime::DAQTicks DAQTime) {
  totRUsageEatTimer_.Start();
  n_hits_in_++;
  //convert to an easier to transport object
  AbsDAQHit h(hashService_->HashFromOMKey(dom), DAQTime);
  
  //turn the cank
  totRUsageHiveTriggerTimer_.Start();
  hiveTrigger_->AddHit(h);
  totRUsageHiveTriggerTimer_.Stop();
  totRUsageEatTimer_.Stop();
//   CollectTriggers();
};

void IceHiveTrigger::AdvanceTime(const hitspooltime::DAQTicks DAQTime) {
  hiveTrigger_->AdvanceTime((hivetrigger::DAQTicks)DAQTime);
  CollectTriggers();
};

void IceHiveTrigger::CollectTriggers() {
  log_debug("CollectTriggers()");
  using namespace hitspooltrigger;
  
  AbsDAQHitSetSequence ht_subEvents = hiveTrigger_->PullSubEvents();
  
  AbsDAQHitSetSequence::const_iterator subEvent = ht_subEvents.begin();
  while (subEvent!= ht_subEvents.end()) {
    if (subEvent->size()>=minEventSize_)
      waitingTriggers_.insert(TriggerWindow(subEvent->begin()->GetDAQTicks(), subEvent->rbegin()->GetDAQTicks()));
    subEvent = ht_subEvents.erase(subEvent);
  }
  
  //merge all overlapping triggers
  if (waitingTriggers_.empty()) {
    return;
  }

  TriggerSet::const_iterator work = waitingTriggers_.begin();
  TriggerSet::const_iterator next = waitingTriggers_.begin(); ++next;
  while (next != waitingTriggers_.end()) {
    if (work->overlaps(*next)) {
      const TriggerWindow tw(work->start, next->end);
      next=waitingTriggers_.erase(next);
      work=waitingTriggers_.erase(work);
      work=waitingTriggers_.insert(work, tw);
    }
    else {
      triggerQueue_.push(*work);
      work = waitingTriggers_.erase(work);
      ++next;
    }
  }
};

void IceHiveTrigger::ReportState() const {};

