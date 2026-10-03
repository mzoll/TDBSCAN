//
// Created by mzoll on 30/08/2026.
//

#ifndef TDBSCAN__COMMON_DEFS_HH
#define TDBSCAN__COMMON_DEFS_HH

#include <format>
#include "common_defs.h"

// ========================= ScalarTime =======================
namespace tdbscan {
bool
ScalarTime_t::operator<(const ScalarTime_t &rhs) const { return value_ < rhs.value_; };

bool
ScalarTime_t::operator==(const ScalarTime_t &rhs) const { return value_ == rhs.value_; };

ScalarTime_t::TimeDiff_t
ScalarTime_t::operator-(const ScalarTime_t &rhs) const { return value_ - rhs.value_; };

//implicit conversion operator for shorthand
ScalarTime_t::operator double() const { return value_; }

//assignment operator
ScalarTime_t &
ScalarTime_t::operator=(const double rhs) {
  value_ = rhs;
  return *this;
};
};  // namespace tdbscan

std::ostream &operator<<(std::ostream &os, const tdbscan::ScalarTime_t &st) {
  return os << std::format("ScalarTime({})", double(st));
};


template<>
struct std::formatter<tdbscan::ScalarTime_t> {
  constexpr auto parse(std::format_parse_context &ctx) { return ctx.begin(); }

  auto format(const tdbscan::ScalarTime_t &st, std::format_context &ctx) const {
    return std::format_to(ctx.out(), "t_{}", static_cast<double>(st));
  };
};


// ========================= Ordinate1d =======================


namespace tdbscan {

template<class t_base, class t_dist>
Position1d<t_base, t_dist>::Distance_t
Position1d<t_base, t_dist>::distance(const Position1d &rhs) const { return rhs.value_ - value_; };

template<class t_base, class t_dist>
Position1d<t_base, t_dist>::Distance_t
Position1d<t_base, t_dist>::magnitude() const { return value_; };

template<class t_base, class t_dist>
Position1d<t_base, t_dist>::Distance_t
Position1d<t_base, t_dist>::abs() const { return this->magnitude(); };

template<class t_base, class t_dist>
bool
Position1d<t_base, t_dist>::operator<(const Position1d &rhs) const { return value_ < rhs.value_; };

template<class t_base, class t_dist>
bool
Position1d<t_base, t_dist>::operator==(const Position1d &rhs) const { return value_ == rhs.value_; };

template<class t_base, class t_dist>
Position1d<t_base, t_dist>::operator double() const { return value_; }

template<class t_base, class t_dist>
Position1d<t_base, t_dist> &
Position1d<t_base, t_dist>::operator=(const double rhs) {
  value_ = rhs;
  return *this;
};

template<class t_base, class t_dist>
std::ostream &operator<<(std::ostream &os, const tdbscan::Position1d<t_base, t_dist> &p1d) {
  return os << std::format("Position1d({})", p1d.getValue());
};
} //namespace tdbscan

template<class t_base, class t_dist>
struct std::formatter<tdbscan::Position1d<t_base, t_dist> > {
  constexpr auto parse(std::format_parse_context &ctx) { return ctx.begin(); }

  auto format(const tdbscan::Position1d<t_base, t_dist> &p1d, std::format_context &ctx) const {
    return std::format_to(ctx.out(), "[x:{}]", p1d.getValue());
  };
};


// --- ContPos1d
namespace tdbscan{
std::ostream &operator<<(std::ostream &os, const tdbscan::ContPos1d &p1d) {
  return os << std::format("ContPos1d({})", p1d.value_);
};
} //namespace tdbscan

template<>
struct std::formatter<tdbscan::ContPos1d> : public std::formatter<tdbscan::Position1d<double, double> > {
//  auto format(const tdbscan::ContPos1d &p1d, std::format_context &ctx) const {
//    return std::formatter<tdbscan::Position1d<double, double> >::format(p1d, ctx);
//  };
};


// ========================= Ordinate3d =======================

namespace tdbscan {

template<class t_base, class t_dist>
Position3d<t_base, t_dist>::Distance_t
Position3d<t_base, t_dist>::distance(const Position3d &rhs) const {
  return sqrt(pow(xord - rhs.xord, 2) + pow(yord - rhs.yord, 2) + pow(zord - rhs.zord, 2));
};

template<class t_base, class t_dist>
Position3d<t_base, t_dist>::Distance_t
Position3d<t_base, t_dist>::magnitude() const { return sqrt(pow(xord, 2) + pow(yord, 2) + pow(zord, 2)); };

template<class t_base, class t_dist>
Position3d<t_base, t_dist>::Distance_t Position3d<t_base, t_dist>::abs() const { return this->magnitude(); };

template<class t_base, class t_dist>
bool Position3d<t_base, t_dist>::operator<(const Position3d &rhs) const {
  return magnitude() < rhs.magnitude() || xord < rhs.xord || yord < rhs.yord || zord < rhs.zord;
};

template<class t_base, class t_dist>
bool Position3d<t_base, t_dist>::operator==(const Position3d &rhs) const {
  return xord == rhs.xord && yord == rhs.yord && zord == rhs.zord;
};

template<class t_base, class t_dist>
std::ostream &operator<<(std::ostream &os, const Position3d<t_base, t_dist> &p3d) {
  return os << std::format("Position3d(x:{}, y:{}, z:{})", p3d.xord, p3d.yord, p3d.zord);
};
}; //namespace tdbscan


template<class t_base, class t_dist>
struct std::formatter<tdbscan::Position3d<t_base, t_dist> > {
  constexpr auto parse(std::format_parse_context &ctx) { return ctx.begin(); }

  auto format(const tdbscan::Position3d<t_base, t_dist> &p3d, std::format_context &ctx) const {
    return std::format_to(ctx.out(), "[x:{}, y:{}, z:{}]", p3d.getXord(), p3d.getYord(), p3d.getZord());
  };
};

// --- ContPos3d
namespace tdbscan{
std::ostream &operator<<(std::ostream &os, const tdbscan::ContPos3d &p3d) {
  return os << std::format("ContPos1d(x:{}, y:{}, z:{})", p3d.xord, p3d.yord, p3d.zord);
};
} //namespace tdbscan

template<>
struct std::formatter<tdbscan::ContPos3d>  : public std::formatter<tdbscan::Position3d<double, double> > {
  //  auto format(const tdbscan::ContPos3d &p3d, std::format_context &ctx) const {
  //    return std::formatter<tdbscan::Position3d<double, double> >::format(p3d, ctx);
  //  };
};


#endif //TDBSCAN__COMMON_DEFS_HH
