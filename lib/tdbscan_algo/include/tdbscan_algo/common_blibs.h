//
// Created by mzoll on 30/08/2026.
//

#ifndef TDBSCAN_COMMON_DEFS_H
#define TDBSCAN_COMMON_DEFS_H


#include <ostream>
#include <format>
#include "tdbscan_algo/absblib.h"
#include "tdbscan_algo/common_defs.h"

namespace tdbscan {
// ========================= BLIB =======================
//make a declaration of the Blib
class Blib3d : public AbsBlib<Position3d, ScalarTime_t> {
public: //type shorthands
  using Ordinate_t = Position3d;
  using Time_t = ScalarTime_t;

protected:
  Ordinate_t pos;
  Time_t time;
public:
  [[nodiscard]] Position3d
  getOrdinate() const
  {return pos;};

  [[nodiscard]] Time_t
  getTime() const
  {return time;};

  [[nodiscard]] Ordinate_t::Distance_t
  distanceTo(const Blib3d& rhs) const
  {return pos.distance(rhs.pos);};

  /// get the time difference
  [[nodiscard]] Time_t::TimeDiff_t
  timeTo(const Blib3d& other) const
  {return other.time - time;};
public: //comparators
  /// define the lesser-operator
  [[nodiscard]] bool
  operator<(const Blib3d& other) const
  { return time < other.time || time == other.time && pos < other.pos; };

  [[nodiscard]] bool
  operator==(const Blib3d& other) const
  { return time == other.time && pos == other.pos; };

  ///constructor
  Blib3d(const Position3d pos, const Time_t time) : pos(pos), time(time) {};

private:
  friend
  std::ostream& operator<< ( std::ostream& os, const Blib3d & b4d);
};


/**
 * A Blib that has a scalar Ordinate
 */
class ScalarBlib : public AbsBlib<Position1d, ScalarTime_t> {
public: //type shorthands
  using Ordinate_t = Position1d;
  using Time_t = ScalarTime_t;

protected:
  Ordinate_t pos;
  Time_t time;
public:
  [[nodiscard]] Position1d
  getOrdinate() const
  {return pos;};

  [[nodiscard]] Time_t
  getTime() const
  {return time;};

  /// get the distance to another blib
  [[nodiscard]] Ordinate_t::Distance_t
  distanceTo(const ScalarBlib& rhs) const
  {return pos.distance(rhs.pos);};

  /// get the time difference to another blib
  [[nodiscard]] Time_t::TimeDiff_t
  timeTo(const ScalarBlib& other) const
  {return other.time - time;};
public: //comparators
  /// define the lesser-operator
  [[nodiscard]] bool
  operator<(const ScalarBlib& other) const
  { return time < other.time || time == other.time && pos < other.pos; };

  [[nodiscard]] bool
  operator==(const ScalarBlib& other) const
  { return time == other.time && pos == other.pos; };

  ///constructor
  ScalarBlib(const Position1d pos, const Time_t time) : pos(pos), time(time) {};

  struct TimeOrder {
    bool operator()(const ScalarBlib& lhs, const ScalarBlib& rhs) const {return double(lhs.getTime()) < double(rhs.getTime());};
  };

public:
  friend
  std::ostream& operator<< ( std::ostream& os, const ScalarBlib & sb);
};

}; //namespace tdbscan

#include "common_blibs.hh"

#endif //TDBSCAN_COMMON_DEFS_H
