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


// ========================= ORDINATE1d =======================

/** make a typedef for an Ordinate that has ONE component
* @tparam t_base the type of the component
* @tparam t_dist the type of the distance measure
*/
template<class t_base, class t_dist>
class Position1d : public Ordinate_t {
public:
  using Distance_t = t_dist;

protected:
  t_base value_{0.};

public:
  /// constructor
  Position1d(const t_base value) : value_(value) {};

  [[nodiscard]] t_base getValue() const {return value_;};

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
  template<class _t_base, class _t_dist>
  friend
  std::ostream &operator<<(std::ostream &os, const Position1d<_t_base, _t_dist> &p1d);
};


/**
* An Ordinate with a single continuous component
*/
class ContPos1d final : public Position1d<double, double> {
  friend
  std::ostream &operator<<(std::ostream &os, const ContPos1d &p1d);
public:
  ContPos1d(const double value) : Position1d(value) {};
};

// =================== Ordinate3d =====================

/** make a typedef for an Ordinate that has THREE component
 *
 * This class is evidently very much similar to a 3vector in math and physics
 * @tparam t_base the type of the components
 * @tparam t_dist the type of the distance measure
 */
template<class t_base, class t_dist>
class Position3d : public Ordinate_t {
public:
  using Distance_t = t_dist;

protected:
  t_base xord, yord, zord;

public:
  /// constructor
  Position3d(const t_base x, const t_base y, const t_base z) : xord(x), yord(y), zord(z) {};

  [[nodiscard]] t_base getXord() const {return xord;}
  [[nodiscard]] t_base getYord() const {return yord;}
  [[nodiscard]] t_base getZord() const {return zord;}

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
  template<class _t_base, class _t_dist>
  friend
  std::ostream &operator<<(std::ostream &os, const Position3d<_t_base, _t_dist> &p3d);
};

/**
* An Ordinate with three continuous components
*/
class ContPos3d final : public Position3d<double, double> {
  friend
  std::ostream &operator<<(std::ostream &os, const ContPos3d &p3d);
public:
  ContPos3d(const double xord, const double yord, const double zord) : Position3d(xord, yord, zord) {};
};
} //namespace tdbscan

#include "common_defs.hh"

#endif //TDBSCAN__COMMON_DEFS_H
