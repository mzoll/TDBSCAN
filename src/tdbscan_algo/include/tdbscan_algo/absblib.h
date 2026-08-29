//
// Created by netsu on 08/08/2026.
//

/**
 * Defines what is an Abstract Blib
 */


#ifndef TDBSCAN_ABSBLIB_H
#define TDBSCAN_ABSBLIB_H


#include <set>
#include <vector>

#include "tdbscan_algo/base_defs.h"

namespace tdbscan {
  /**
   * copies needs to be lightweight.
   * @tparam tOrdinate_t an ordinate class, aka a position
   * @tparam tTime_t a timelike class
   */
  template <class tOrdinate_t, class tTime_t>  //tOrdinate adheres to Position_t, tTime_t aheres to Time_t
  class AbsBlib {
  public:
    using Ordinate_t = tOrdinate_t;
    using Time_t = tTime_t;
  public:  //methods
    [[nodiscard]]
    Ordinate_t
    getOrdinate() const;

    [[nodiscard]]
    Time_t
    getTime() const;
  public: //comparators
    [[nodiscard]]
    typename Ordinate_t::Distance_t
    distanceTo(const AbsBlib& other) const;

    /// get the time difference to rhs
    [[nodiscard]]
    typename Time_t::TimeDiff_t
    timeDiff(const AbsBlib& other) const;

    [[nodiscard]]
    bool
    operator<(const AbsBlib& other) const;

    [[nodiscard]]
    bool
    operator==(const AbsBlib& other) const;
  };

  template <typename tPosition>
  using AbsBlibSet = std::set< AbsBlib<tPosition, Time_t> >;
}

#endif //TDBSCAN_ABSBLIB_H

