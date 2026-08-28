//
// Created by netsu on 07/08/2026.
//

#ifndef TDBSCAN_TDBSCAN_ALGO_H
#define TDBSCAN_TDBSCAN_ALGO_H

#include <vector>

#include "base_defs.h"
#include "tdbcluster.h"

namespace tdbscan {

namespace detail {

template <class tBlib>
bool CausallyConnected(
  const tBlib& h1,
  const tBlib& h2,
  const Connector<tBlib>& connector)
{
  if (h1.GetTime() > h2.GetTime())
    return CausallyConnected(h2, h1, connector); //recursive call to enforce time-order at this point
  return connector.Connected(h1, h2);
}
} // namespace tdbscan::detail

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
public: //typedefs: some internal definitions and shorthands
  //  SET_LOGGER("HiveSplitter");

  /// A set of hits, time-order is enforced automatically
  typedef std::set<tBlib> BlibSet;

  struct BlibSetTimeOrder {
    /// implement the order principle
    bool operator()(const BlibSet &lhs, const BlibSet &rhs) const;
  };

  ///a time-ordered sequence of time-ordered HitSets
  typedef std::set<BlibSet, BlibSetTimeOrder> BlibSetSequence;

  typedef Connector<tBlib> Connector_t;
  typedef CausalCluster<tBlib> CausalCluster_t;

private: // internal state
  //==================
  // Properties
  //==================

  /// the time of the algo
  Time_t sync_time;

  /// all in-progress causal clusters
  std::list<CausalCluster_t> active_clusters_;
  /// all concluded clusters
  std::list<CausalCluster_t> clusters_;

private: //parameters
  //========================
  // Configurable Parameters
  //========================
  /// PARAM: A parameter-set to run on
  TDBScan_ParameterSet params_;
  /// PARAM: this defines the 'physics' at play
  Connector_t* connector_;

public: //interface
  /**
   * Constructor from a ParameterSet and Connector
   * @param params Collection of Parameters which govern the algorithm
   * @param connector Pointer to a Connector, which facilítates the comparison of hits
   */
  TDBScan_Algo(
    const TDBScan_ParameterSet& params,
    const Connector_t* connector);

  /** @brief ACTION
   * Perform the Splitting feeding it a series of Hits
   * Than call the Main-Algorithm and iterate over all the hits, try to form clusters of Hits and add Hits to them, discard them or form a new cluster.
   * When Hit-Series is exhausted and Clusters/Subevents have been found write them to the data-stream (either as separate frames or into the same frame).
   * Save the SplitCount also if that is wanted,
   * Push everything back into te pipeline.
   * Clean-up the memory.
   * @tparam tBlibContainer any type of iterable object containing tBlibs
   * @param inhits the hits to process
   * @return a series of hits, which are the subevents (time-order in sequence and in hit-order)
   */
  template <class tBlibContainer>
  BlibSetSequence Process (const tBlibContainer& blibs);

  BlibSetSequence Process (const BlibSet& blibs);

public: // probe of internal state
  /// Get the time until which the result is static and no active hits are percolating in the algorithm/clusters
  [[nodiscard]]
  Time_t FinalizedUntil() const;

  //probe into the algorithm
  bool CausallyConnected(const tBlib& b1, const tBlib& b2) const;



protected: // --- THE REAL MACHINERY ---
  //===================
  // Internal Methods
  //===================

  /**The main driver for the entire algorithm:
   * Adds a new hit to all clusters with which it is connected (including subsets of existing clusters).
   * By 'advancing' the clusters this function also causes subevents to be built when possible.
   * @param b the blib to add
   */
  void NextBlib(const tBlib &b);

private: //work on *clusters*



  /// advance the cluster forward in time TODO more specific
  void AdvanceInTime(tdbscan::CausalCluster<tBlib>& c, Time_t timep);

  /** Attempt to add Hit h to existing cluster c, or to the subset of c with which it is connected by enough hits in c
   * to meet the multiplicity condition.
   * @param c the cluster to add to
   * @param b the blib to add
   */
  void InsertHit(
    tdbscan::CausalCluster<tBlib>& c,
    const tBlib& b);  // This needs modification

private: //things that work on the surface of Clusters, but do not change the internal state
  /// Does Hit h connect to CausalCluster
  [[nodiscard]]
  bool connectsTo(const tBlib& b, const CausalCluster<tBlib>& c1) const;

  void evalEstablished(CausalCluster<tBlib> c) const {
    if (c.count() >= params_.multiplicity)
      c.established = true;
  }

  /** get the CausalCluster of all active hits within this cluster which can be considered connected
   * \param h the Hit to check against
   */
  [[nodiscard]]
  CausalCluster<tBlib> getConnectedSubCluster(
    const tBlib &h, const CausalCluster<tBlib>& c1) const;

  ///insert an hit and advance the cluster
  void insertActiveHit(const tBlib &h, CausalCluster<tBlib>& c) {
    c.insertActiveHit(h);
    if (c.active_doms.size()>=params_.multiplicity) {
      c.established=true;
    }
  }

  void Finalize();

};
};

#include "tdbscan_algo.hh"

#endif //TDBSCAN_TDBSCAN_ALGO_H
