//
// Created by netsu on 07/08/2026.
//

#ifndef TDBSCAN_TDBCLUSTER_H
#define TDBSCAN_TDBCLUSTER_H

#include <limits>
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

namespace detail {


  ///sufficient overlap in set1 and set2 by hits on 'multiplicity' many DOMs with 'multiplicityTimeWindow'
  //(agnostic)
  template <class tPosition>
  bool CausallyOverlaps (
    const AbsBlibSet<tPosition>& set1,
    const AbsBlibSet<tPosition>& set2,
    const unsigned int multiplicity,
    const Time_t multiplicityTimeWindow);
} //namespace detail

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

  public: //properties
    ///The ordered set of hits within this cluster
    BlibSet blibs_;

    enum Status {
      EMPTY = 0,
      EMERGING = 10,
      GROWING = 20,
      DYING = 30,
      CONCLUDED = 40,
    } status = EMPTY;


  public:
    // Constructor
    CausalCluster();
    ///adhoc constructor from single blib
    CausalCluster(const tBlib &h );
    ///adhoc constructor series of hits
    CausalCluster(const std::set<tBlib> &bset);
    /// this should be made in a proper copy constructor
    CausalCluster(const CausalCluster& cc);

  protected: //methods (altering)
    ///Add a new hit to the cluster
    ///\param h The hit to add
    void insertBlib(const tBlib &h);
    /// Take all hits from the cluster and add them to its own
    /// @param c the cluster to be merged
    void copyHits(const CausalCluster& c);

  public: //methods (inert)
    /// get hits of this cluster
    [[nodiscard]] const BlibSet& getHits() const;
    ///Finds the time of the earliest hit in this cluster
    /// @return The earliest hit time or minus infinity if the cluster is empty
    [[nodiscard]] typename tBlib::Time_t
    getEarliestTime() const;
    ///Finds the time of the latest hit in this cluster
    /// @return The latest hit time or infinity if the cluster is empty
    [[nodiscard]] typename tBlib::Time_t
    getLatestTime() const;

    [[nodiscard]] uint64_t
    nHitsWithinTimeWindow(const typename tBlib::Time_t earliest, const typename tBlib::Time_t latest) const;


    ///is this cluster established
    [[nodiscard]] bool isEstablished() const;

    /// ist this cluster concluded (no hits can be added)
    [[nodiscard]] bool isConcluded() const;

    ///Test whether the hits in this are a subset those of super
    /// @param super cluster with a series of hits which might be a superset
    /// @return true, if this is a subset of super
    [[nodiscard]] bool isSubsetOf(const CausalCluster& super) const;

    ///Test whether the hits in this are a subset those of super
    /// @param sub cluster with a series of hits which might be a superset
    /// @return true, if this is a subset of super
    [[nodiscard]] bool isSupersetOf(const CausalCluster& sub) const;

    /// Test whether two Clusters contain the same hits
    [[nodiscard]] bool isConcruent(const CausalCluster& c2) const;

    /// The number of blibs in both Clusters
    [[nodiscard]] unsigned int nOverlap(const CausalCluster& c2) const;




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


}// namespace tdbscan

#include "tdbcluster.hh"

#endif //TDBSCAN_TDBCLUSTER_H
