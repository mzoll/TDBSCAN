//
// Created by mzoll on 30/08/2026.
//

#ifndef TDBSCAN__COMMON_DEFS_H
#define TDBSCAN__COMMON_DEFS_H

#include <cmath>
#include <format>
#include "tdbscan_algo/base_defs.h"


namespace tdbscan {
// ========================= TIME =======================

/** make a typedef for what is the notion of Time;
* On first principles time is a continuous monotonic increasing variable
*/
class ScalarTime_t final : Time_t {
public:
  typedef double TimeDiff_t;

private:
  double value_{0.};

public:
  /// default constructor: for convenience
  ScalarTime_t() : value_(0.) {};
  ScalarTime_t(const double value) : value_(value) {};

  inline
  bool
  operator<(const ScalarTime_t &rhs) const;

  inline
  bool
  operator==(const ScalarTime_t &rhs) const;

  inline
  TimeDiff_t
  operator-(const ScalarTime_t &rhs) const;

  //implicit conversion operator for shorthand
  operator double() const;

  //assignment operator
  ScalarTime_t &operator=(const double rhs);

  static constexpr double min() { return -std::numeric_limits<double>::infinity(); };
  static constexpr double max() { return std::numeric_limits<double>::infinity(); };

private:
  friend
  std::ostream &operator<<(std::ostream &os, const ScalarTime_t &st);
};


// ========================= ORDINATE =======================

/** make a typedef for what is the notion of Time;
* On first principles time is a continuous monotonic increasing variable
*/
class Position1d final : public tdbscan::Ordinate_t {
public:
  typedef double Distance_t;

private:
  double value_{0.};

public:
  /// constructor
  Position1d(const double value) : value_(value) {};

  /// the distance to another position
  [[nodiscard]]
  Distance_t
  distance(const Position1d &rhs) const;

  [[nodiscard]]
  Distance_t magnitude() const;;

  [[nodiscard]]
  Distance_t abs() const;

  [[nodiscard]]
  bool operator<(const Position1d &rhs) const;

  [[nodiscard]]
  bool operator==(const Position1d &rhs) const;

public: //convenience
  ///implicit conversion operator for shorthand
  explicit operator double() const;

  /// assignment operator
  Position1d &operator=(double rhs);

private:
  friend
  std::ostream &operator<<(std::ostream &os, const Position1d &p1d);
};


/// make a definition of a Point in 3d space
class Position3d final : public tdbscan::Ordinate_t {
public:
  typedef double Distance_t;

public:
  double xord, yord, zord;

public:
  /// constructor
  Position3d(const double x, const double y, const double z) : xord(x), yord(y), zord(z) {};

  /// get the distance with a partner object
  [[nodiscard]]
  Distance_t
  distance(const Position3d &rhs) const;

  [[nodiscard]]
  Distance_t magnitude() const;

  [[nodiscard]]
  Distance_t abs() const;

  [[nodiscard]]
  bool operator<(const Position3d &rhs) const;

  [[nodiscard]]
  bool operator==(const Position3d &rhs) const;

private:
  friend
  std::ostream &operator<<(std::ostream &os, const Position3d &p3d);
};
} //namespace tdbscan

#include "common_defs.hh"

#endif //TDBSCAN__COMMON_DEFS_H
