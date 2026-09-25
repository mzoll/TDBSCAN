//
// Created by mzoll on 07/08/2026.
//

#ifndef TDBSCAN_TDBCLUSTER_H
#define TDBSCAN_TDBCLUSTER_H

#include <list>
#include <set>
#include <cstdint>

#include "absblib.h"
#include "connector.h"



// --- MACHINERY PARTS ---
//======================================
// Data Structures and Helper Functions
//======================================
namespace tdbscan {
  template<typename tBlib>
  using BlibSet = std::set<tBlib>;

  template<typename tBlib>
  using BlibList = std::list<tBlib>;

  /** An object which keeps track of a group of hits which are (mostly) causally connected to each other,
   * and the number of distinct DOMs on which those hits occurred.
   * There is some internal bookkeeping and intelligence on that
   *
   * @tparam tBlib the concrete Blib class used; needs to derive from AbsBlib
   */
  template <class tBlib>
  class CausalCluster{
  public:
    //friend class TDBScan_Algo;  this is a forward declaration for friend access

    using BlibSet = std::set<tBlib>;
    using BlibList = std::list<tBlib>;

    using tTime = tBlib::Time_t;

  public: //properties
    ///The ordered set of hits within this cluster
    BlibSet blibs_;

  public:
    // Constructor
    CausalCluster();
    ///adhoc constructor from single blib
    CausalCluster(const tBlib &h);
    ///adhoc constructor series of hits
    CausalCluster(const std::set<tBlib> &bset);
    /// this should be made in a proper copy constructor
    CausalCluster(const CausalCluster& cc);

  public: //methods (altering)
    ///Add a new hit to the cluster
    /// @param h The hit to add
    void insertBlib(const tBlib &h);
    /// Take all hits from the cluster and add them to its own
    /// @param c the cluster to be merged
    void copyBlibs(const CausalCluster& c);

  public: //methods (inert)
    /// get hits of this cluster
    [[nodiscard]] const BlibSet& getHits() const;
    ///Finds the time of the earliest hit in this cluster
    /// @return The earliest hit time or minus infinity if the cluster is empty
    [[nodiscard]] tTime
    getEarliestTime() const;
    ///Finds the time of the latest hit in this cluster
    /// @return The latest hit time or infinity if the cluster is empty
    [[nodiscard]] tTime
    getLatestTime() const;

    [[nodiscard]] uint64_t
    nHitsWithinTimeWindow(tTime earliest = tTime::min(), tTime latest = tTime::max()) const;

    ///Test whether the hits in this are a subset those of super
    /// @param super cluster with a series of hits which might be a superset
    /// @return true, if this is a subset of super
    [[nodiscard]] bool isSubsetOf(const CausalCluster& super) const;

    ///Test whether the hits in this are a subset those of super
    /// @param sub cluster with a series of hits which might be a superset
    /// @return true, if this is a subset of super
    [[nodiscard]] bool isSupersetOf(const CausalCluster& other) const;

    /// Test whether two Clusters contain the same hits
    [[nodiscard]] bool isConcruent(const CausalCluster& other) const;

    /// The number of blibs in both Clusters
    [[nodiscard]] unsigned int nOverlap(const CausalCluster& other) const;

    [[nodiscard]]
    inline
    bool empty() const;

    [[nodiscard]]
    inline
    uint64_t count() const;
  };


  /// an alias shorthand
  template<typename T>
  using CausalClusterList = std::list<CausalCluster<T> >;

} // namespace tdbscan

#include "tdbcluster.hh"

#endif //TDBSCAN_TDBCLUSTER_H
