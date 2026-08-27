/**
 * \file HiveTrigger.h
 *
 * (c) 2012 the IceCube Collaboration
 *
 * $Id: HiveTrigger.h 99900 2013-02-26 10:10:43Z mzoll $
 * \version $Revision: 99900 $
 * \date $Date: 2013-02-26 11:10:43 +0100 (Tue, 26 Feb 2013) $
 * \author Marcel Zoll <marcel.zoll@fysik.su.se>
 *
 * The central algorithm to split ResponseSeriesMaps by arguments of clustering
 */

#ifndef HIVETRIGGER_H
#define HIVETRIGGER_H

#include <limits>
#include <list>
#include <map>
#include <set>
#include <deque>

#include "ToolZ/OMKeyHash.h"
#include "ToolZ/HitSorting.h"
#include "IceHiveZ/internals/Connector.h"

namespace hivetrigger {
  
  ///make a typedef for what is the notion of Time in I3Units::Ns
  typedef double Time;
  ///make a typedef for the internal state of affairs: DAQTicks
  typedef int64_t DAQTicks;
  ///convert DAQTicks to I3UnitsNs
  DAQTicks NsToTicks(const Time ns);
  ///convert I3UnitsNs to DAQTicks
  Time TicksToNs(const DAQTicks ticks);
  
  /// A set of parameters that steer HiveTrigger
  struct HiveTrigger_ParameterSet{
    /// PARAM: Required multiplicity of connected !DOMs! with any hit within the time-window for to be accected to the cluster
    size_t multiplicity;
    /// PARAM: Time span within which the multiplicity requirement must be met
    Time multiplicityTimeWindow;
    /// PARAM: Connect all hits on same DOM up to this time limit after the initial hit regardlessly; deactivate by NAN
    Time acceptTimeWindow;
    /// PARAM: Reject all hits on same DOM from to this time limit after the initial hit regardlessly; deactivate by INF
    Time rejectTimeWindow;
    /// PARAM: pointer to the ConnectorBlock providing DOM to DOM and Hit connections
    ConnectorBlockPtr connectorBlock;
    /// PARAM: number of overlapping !DOMs! required for (partial)subevents to be merged
    size_t mergeOverlap;
    
    ///constructor
    HiveTrigger_ParameterSet();
  };

// --- MACHINERY PARTS ---
//======================================
// Data Structures and Helper Functions
//======================================
namespace detail {
  
  ///enforces timeorder in hits, which is important [h1 should be earlier than h2]
  bool CausallyConnected(
    const AbsDAQHit& h1,
    const AbsDAQHit& h2,
    const ConnectorBlockConstPtr& connectorBlock);
  
  ///sufficent overlap in set1 and set2 by hits on 'multiplicity' many DOMs with 'multiplicityTimeWindow'
  bool CausallyOverlaps (
    const AbsDAQHitSet& set1,
    const AbsDAQHitSet& set2,
    const size_t multiplicity,
    const DAQTicks multiplicityTimeWindow);
  
  ///An object which keeps track of a group of hits which are (mostly) causally connected to each other,
  ///and the number of distinct DOMs on which those hits occurred
  class CausalCluster{
  private: //param
    /// a major steering set of parameters
    const HiveTrigger_ParameterSet* params;
    
  private: //properties
    ///the latest time to which this cluster is syncronized
    DAQTicks sync_time;
    ///The ordered queue of hits within this cluster which are still within the time window of the current time
    AbsDAQHitSet active_hits;
    ///Keeps track of the number of hits on each of the doms present in this cluster, keys are dom indices
    typedef std::map<CompactHash, size_t> DOMHitCount;
    DOMHitCount active_doms;
    ///the DOMs and time of their first hit
    typedef std::map<CompactHash, DAQTicks> DOMHitTimes;
    DOMHitTimes firstHitTimes;
    ///The hits which have formed a group surpassing the multiplicity and are now outside the time window
    AbsDAQHitSet concluded_hits;
    ///Whether the multiplicity condition is met
    bool established;

  public://methods
    ///constructor
    ///\param p the parameter set, which contains essential information when to connect hits
    CausalCluster(const HiveTrigger_ParameterSet* p);
    ///The active hits of this cluster has enough overlap with this hit, so that it should be considered connected
    ///\param h The hit to check
    bool connectsTo(const AbsDAQHit &h) const;
    ///Add a new hit to the cluster
    ///\param h The hit to add
    void insertActiveHit(const AbsDAQHit &h);
    ///Take all hits in other's concluded_hits list and merge them into this cluster's concluded_hits list
    ///\param c the cluster to be merged
    void takeConcludedHits(const CausalCluster& c);
    ///get the concluded hits of this cluster
    const AbsDAQHitSet& getActiveHits() const;
    ///get the concluded hits of this cluster
    const AbsDAQHitSet& getConcludedHits() const;
    ///get the concluded hits of this cluster
    const DOMHitTimes& getFirstHitTimes() const;
    ///get the latest hit, i.e. the most recently added hit
    const AbsDAQHit& getLatestActiveHit() const;
    ///get the CausalClaster of all active hits within this cluster which can be considered connected
    ///\param h the Hit to check against
    CausalCluster getSubCluster(const AbsDAQHit &h) const;
    ///Move this cluster forward in time to t, dropping hits which are no longer within the time window,
    ///the request to merge clusters is accounted for
    ///\param time The current time to which the cluster should be moved
    void advanceInTime(const DAQTicks time);
    ///Finds the time of the earliest hit in this cluster
    ///\return The earliest hit time or infinity if the cluster is empty
    DAQTicks getEarliestTime() const;
    ///Finds the time of the latest hit in this cluster
    ///\return The latest hit time or infinity if the cluster is empty
    DAQTicks getLatestTime() const;
    ///is this cluster still active; thus can there still be found connected hits?
    bool isActive() const;
    ///is this cluster still active; thus can there still be found connected hits?
    bool isEstablished() const;
    ///Test whether the hits in sub are a subset of super
    ///\param super cluster with a series of hits which might be a superset
    ///\return true, if sub is a subset of super
    bool isSubsetOf(const CausalCluster& super) const;
  };

  typedef std::list<CausalCluster> CausalClusterList;
  
}// namespace detail
}// namespace hivetrigger


///The main splitter module
class HiveTrigger {
  SET_LOGGER("HiveTrigger");
private: //properties/internals
  //==================
  // Properties / internals
  //==================
  ///all in-progress causal clusters
  hivetrigger::detail::CausalClusterList clusters_;
  ///temporary storage for causal clusters generated while adding a single hit
  hivetrigger::detail::CausalClusterList newClusters_;
  ///all in-progress subevents
  AbsDAQHitSetList partialSubEvents_;
public: //exposed internals
  ///set of completed subevents which are in time-order (in every aspect)
  AbsDAQHitSetSequence subEvents_;
protected: //parameters
  //========================
  // Configurable Parameters
  //========================
  /// PARAM: A parameter-set to run on
  hivetrigger::HiveTrigger_ParameterSet params_;
  
public: //interface
  //===================
  // Interface
  //===================
  /// Constructor from a ParameterSet
  HiveTrigger(const hivetrigger::HiveTrigger_ParameterSet& params);
  
  /// Get the time until which the result is static and no active hits are perculating in the algorithm/clusters
  hivetrigger::DAQTicks FinalizedUntil() const;
//   hivetrigger::Time FinalizedUntil() const;
  
  /// retrieve all finished SubEvents
  AbsDAQHitSetSequence PullSubEvents();

  /**The main driver for the entire algorithm:
   * Adds a new hit to all clusters with which it is connected (including subsets of existing clusters).
   * By 'advancing' the clusters this function also causes subevents to be built when possible.
   * @param h the hit to add
   */
  void AddHit(const AbsDAQHit &h);

private:
  //================
  // Internal Methods
  //================
  
  /** Attempt to add Hit h to existing cluster c, or to the subset of c with which it is connected by enough hits in c
   * to meet the multiplicity condition.
   * @param c the cluster to add to
   * @param h the hit to add
   * @return true, if h was added to c, or to a new subset of c;
   *	       false, if h was not placed in any cluster
   */
  bool AddHitToCluster(hivetrigger::detail::CausalCluster& c,
                       const AbsDAQHit& h);

  /** Inserts a cluster of hits into the set of subevents, after merging it with any existing subevents
   * with which it shares at least one hit
   * This function also moves any subevents which can no longer grow into the finished subevent collection.
   * @param newSet the set to add
   */
  void AddSubEvent(AbsDAQHitSet& newSet);

public:
  /** Push all Subevents which are finalized into the concluded Subevents
   * @param earliestTime the time until which subevents can be pushed
   */
  void PushEvents(const hivetrigger::DAQTicks earliestTick);
  void PushEvents(const hivetrigger::Time earliestTime);

  /// advance the existing clusters to this time; no hits can be retrocatively inserted before this time
  /// \param time in DAQ ticks
  void AdvanceTime(const hivetrigger::DAQTicks tick);
  void AdvanceTime(const hivetrigger::Time time);

  /** Pushes all hits through the clusters and completes all subevents,
   * on the assumption that no more future hits will be added.
   */
  void FinalizeSubEvents();
};


//===========================================
//============== IMPLEMENTATION =============
//===========================================


inline
hivetrigger::DAQTicks 
hivetrigger::NsToTicks(const hivetrigger::Time ns) 
  {return (hivetrigger::DAQTicks)(ns*10.);};

inline
hivetrigger::Time
hivetrigger::TicksToNs(const hivetrigger::DAQTicks ticks) 
  {return (hivetrigger::Time)(ticks/10.);};

// inline
// hivetrigger::Time HiveTrigger::FinalizedUntil() const
//   {return TicksToNs(FinalizedUntil());};

inline
void HiveTrigger::PushEvents(const hivetrigger::Time earliestTime)
  {return PushEvents(hivetrigger::NsToTicks(earliestTime));};
  
inline
void HiveTrigger::AdvanceTime(const hivetrigger::Time time)
  {return PushEvents(hivetrigger::NsToTicks(time));}

#endif
