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
   * @tparam tPosition
   */
  template <class tPosition, class tTime_t>  //tPosition adheres to Position_t
  class AbsBlib {
  public:  //methods
    [[nodiscard]]
    tPosition GetPosition() const;

    [[nodiscard]] tTime_t GetTime() const;
  public: //comparators
    [[nodiscard]] Distance_t GetDistance(const AbsBlib& other) const;

    /// get the time difference to rhs
    [[nodiscard]] Timediff_t TimeDiff(const AbsBlib& other) const;

    [[nodiscard]] bool operator<(const AbsBlib& other) const;
  };

  template <typename tPosition>
  using AbsBlibSet = std::set< AbsBlib<tPosition, Time_t> >;
}

#endif //TDBSCAN_ABSBLIB_H

