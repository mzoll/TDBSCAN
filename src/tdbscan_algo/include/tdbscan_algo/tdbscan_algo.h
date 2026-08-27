//
// Created by netsu on 07/08/2026.
//

#ifndef TDBSCAN_TDBSCAN_ALGO_H
#define TDBSCAN_TDBSCAN_ALGO_H

#include <vector>

#include "base_defs.h"
#include "tdbcluster.h"

namespace tdbscan {
  /// A set of parameters that steer HiveSplitter
  struct TDBScan_ParameterSet{
    /// PARAM: Required multiplicity of connected !DOMs! with any hit within the time-window for to be accepted to the cluster
    unsigned int multiplicity;
    /// PARAM: Time span within which the multiplicity requirement must be met
    Time_t multiplicityTimeWindow;
    /// PARAM: Connect all hits on same DOM up to this time limit after the initial hit regardlessly; deactivate by NAN
    Time_t acceptTimeWindow;
    /// PARAM: Reject all hits on same DOM from to this time limit after the initial hit regardlessly; deactivate by INF
    Time_t rejectTimeWindow;
    /// PARAM: number of overlapping !DOMs! required for (partial)subevents to be merged into a super set
    unsigned int mergeOverlap;

    ///constructor
    TDBScan_ParameterSet();
  };


  /**
   * The main algorithm class
   * Needs to be configured with a connector and a parameter set then can be fed with a sequence of hits
   */
  template <class tBlib>
  class TDBScan_Algo {
    //  SET_LOGGER("HiveSplitter");

    typedef std::set<tBlib> BlibSet;

    struct BlibSetTimeOrder {
      /// implement the order principle
      bool operator()(const BlibSet &lhs, const BlibSet &rhs) const;
    };

    //a time-ordered sequence of time-ordered HitSets
    typedef std::set<BlibSet, BlibSetTimeOrder> BlibSetSequence;


  private: // internal state
    //==================
    // Properties
    //==================

    /// the time of the algo
    Time_t sync_time;

    /// all in-progress causal clusters
    CausalClusterList<tBlib> active_Clusters_;
    /// all concluded clusters
    CausalClusterList<tBlib> clusters_;

  private: //parameters
    //========================
    // Configurable Parameters
    //========================
    /// PARAM: A parameter-set to run on
    TDBScan_ParameterSet params_;
    /// PARAM: this defines the 'physics' at play
    ConnectorPtr connector_;

  public: //interface
    /**
     * Constructor from a ParameterSet and Connector
     * @param params Collection of Parameters which govern the algorithm
     * @param connector Pointer to a Connector, which facilítates the comparison of hits
     */
    TDBScan_Algo(
      const TDBScan_ParameterSet& params,
      const ConnectorPtr& connector);

    /** @brief ACTION
     * Perform the Splitting feeding it a series of Hits
     * Than call the Main-Algorithm and iterate over all the hits, try to form clusters of Hits and add Hits to them, discard them or form a new cluster.
     * When Hit-Series is exhausted and Clusters/Subevents have been found write them to the data-stream (either as separate frames or into the same frame).
     * Save the SplitCount also if that is wanted,
     * Push everything back into te pipeline.
     * Clean-up the memory.
     * @tparam tBlibContainer any type of iteratable object containing tBlibs
     * @param inhits the hits to process
     * @return a series of hits, which are the subevents (time-order in sequence and in hit-order)
     */
    template <class tBlibContainer>
    BlibSetSequence Process (const tBlibContainer& inhits);

  public: // probe of internal state
    /// Get the time until which the result is static and no active hits are percolating in the algorithm/clusters
    [[nodiscard]]
    Time_t FinalizedUntil() const;

  private: // --- THE REAL MACHINERY ---
    //===================
    // Internal Methods
    //===================

    /**The main driver for the entire algorithm:
     * Adds a new hit to all clusters with which it is connected (including subsets of existing clusters).
     * By 'advancing' the clusters this function also causes subevents to be built when possible.
     * @param h the hit to add
     */
    void AddHit(const tBlib &h);


  private: //work on *clusters*
    void AdvanceInTime(tdbscan::CausalCluster<tBlib>& c, Time_t timep);

    /** Attempt to add Hit h to existing cluster c, or to the subset of c with which it is connected by enough hits in c
     * to meet the multiplicity condition.
     * @param c the cluster to add to
     * @param h the hit to add
     * @return true, if h was added to c, or to a new subset of c;
     *	       false, if h was not placed in any cluster
     */
    bool AddHitToCluster(
      tdbscan::CausalCluster<tBlib>& c,
      const tBlib& h);  // This needs modification

  private: //things that work on the surface of Clusters, but do not change the internal state
    /// Does Hit h connect to CausalCluster
    [[nodiscard]]
    bool connectsTo(const tBlib& h, const CausalCluster<tBlib>& c1) const;

    void evalEstablished(CausalCluster<tBlib> c) const {
      if (c.count() >= params_.multiplicity)
        c.established = true;
    }

    /** get the CausalCluster of all active hits within this cluster which can be considered connected
     * \param h the Hit to check against
     */
    [[nodiscard]]
    CausalCluster<tBlib> getConnectedSubCluster(const tBlib &h, const CausalCluster<tBlib>& c1) const;

    ///insert an hit and advance the cluster
    void insertActiveHit(const tBlib &h, CausalCluster<tBlib>& c) {
      c.insertActiveHit(h);
      if (c.active_doms.size()>=params_->multiplicity) {
        c.established=true;
      }
    }

    ///check if this cluster is still active
    bool isActive(const CausalCluster<tBlib>& c) const {
      if (c.hasActiveHits())
        return true;

      if (params_.acceptTimeWindow <= params_.multiplicityTimeWindow) {
        //then only active hits can accept more hits
        return false;
      }
      //need to look into the time of every first hit on any each DOM if further hits can be accepted
      const auto latest_accept_time = c.sync_time - params_.acceptTimeWindow;
      for (const auto& fht, firstHitTimes) {
        if (fht.second > latest_accept_time)
          return true;
      }
      return false;
    }

    ///Move this cluster forward in time to t, dropping hits which are no longer within the time window,
    ///\param time The current time to which the cluster should be moved
    void advanceInTime(const CausalCluster& c, const Time_t time) {
      while (c.hasActiveHits()) {
        const auto h = c.active_hits.begin();
        if (time > h->GetTime()+ params_.multiplicityTimeWindow) {
          //the hit is no longer active, thus
          //decrement the number of hits on the DOM where h occurred
          if ((--c.active_doms[h->GetDOMIndex()])<=0) //NOTE TODO do we need to bother with this after the cluster is established, and this is probably not checked anymore?
            c.active_doms.erase(h->GetDOMIndex());

          //if the multiplicity threshold was met include h in the finished cluster
          if (c.established) {
            //insert the hit
            c.concluded_hits.insert(c.concluded_hits.end(),*h);
          }
          else { //hit is about to be discarded
            //sync up the firsthit-time map
            if (!c.active_doms.count(h->GetDOMIndex()))
              c.firstHitTimes.erase(h->GetDOMIndex());
            else {
              //check for the next hit on the same DOM which is still active and take its time instead
              BOOST_FOREACH(const AbsHit& hh, c.active_hits) {
                if (h->GetDOMIndex() == hh.GetDOMIndex())
                  c.firstHitTimes[h->GetDOMIndex()] = hh.GetTime();
              }
              //NOTE by this shift some inconsitency is introduced of the connections between hits in the cluster
              //however the merging of Clusters in the HiveSplitter will bring this all in sync again
            }
          }
          c.active_hits.erase(h);
        }
        else
          break;
      }
      //the cluster is now synced to this time
      c.sync_time=time;
    }

  };

  ///specialization for already (time-)sorted hits
  template <>
  AbsHitSetSequence HiveSplitter::Split (const AbsHitSet& inhits);



}

#include "tdbscan_algo.hh"

#endif //TDBSCAN_TDBSCAN_ALGO_H
