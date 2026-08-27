//
// Created by netsu on 08/08/2026.
//

/**
 * Defines what is an Abstract Blib
 */


#ifndef TDBSCAN_ABSBLIB_H
#define TDBSCAN_ABSBLIB_H


#include <vector>

#include "base_defs.h"

namespace tdbscan {
  /**
   * copies needs to be lightweight.
   * @tparam tPosition
   */
  template <class tPosition>  //tPosition adheres to Position_t
  class AbsBlib {
  public:  //methods
    [[nodiscard]] virtual tPosition GetPosition() const = 0;

    [[nodiscard]] virtual Time_t GetTime() const = 0;
  public: //comparators
    [[nodiscard]] virtual Distance_t GetDistance(const AbsBlib& rhs) const = 0;

    /// get the time difference to rhs
    [[nodiscard]] virtual Time_t TimeDiff(const AbsBlib& rhs) const = 0;
  };

  template <typename tPosition>
  using AbsBlibSet = std::set< AbsBlib<tPosition> >;
}

#endif //TDBSCAN_ABSBLIB_H

