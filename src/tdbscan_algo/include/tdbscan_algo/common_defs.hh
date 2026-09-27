//
// Created by mzoll on 30/08/2026.
//

#ifndef TDBSCAN__COMMON_DEFS_HH
#define TDBSCAN__COMMON_DEFS_HH

#include <format>
#include "common_defs.h"

// ========================= TIME =======================
namespace tdbscan {

bool
ScalarTime_t::operator<(const ScalarTime_t &rhs) const
{return value_ < rhs.value_;};

bool
ScalarTime_t::operator==(const ScalarTime_t &rhs) const
{return value_ == rhs.value_;};

ScalarTime_t::TimeDiff_t
ScalarTime_t::operator-(const ScalarTime_t &rhs) const
{return value_ - rhs.value_;};

//implicit conversion operator for shorthand
ScalarTime_t::operator double() const { return value_; }

//assignment operator
ScalarTime_t&
ScalarTime_t::operator=(const double rhs)
{value_=rhs; return *this; };
};

std::ostream& operator<< ( std::ostream& os, const tdbscan::ScalarTime_t & st) {
	return os << std::format("ScalarTime({})", double(st));
};


template <>
struct std::formatter<tdbscan::ScalarTime_t> {
  constexpr auto parse(std::format_parse_context & ctx) {return ctx.begin();}

  auto format(const tdbscan::ScalarTime_t& st, std::format_context &ctx) const {
    return std::format_to(ctx.out(), "t_{}", static_cast<double>(st));
  };
};


namespace tdbscan {

/// the distance to another position
Position1d::Distance_t
Position1d::distance(const Position1d& rhs) const
{ return rhs.value_ - value_;};

Position1d::Distance_t
Position1d::magnitude() const
{return value_;};

Position1d::Distance_t
Position1d::abs() const
{return this->magnitude();};

bool
Position1d::operator<(const Position1d &rhs) const
{return value_ < rhs.value_;};

bool
Position1d::operator==(const Position1d &rhs) const
{return value_ == rhs.value_;};


Position1d::operator double() const { return value_; }

Position1d&
Position1d::operator=(const double rhs)
	{value_=rhs; return *this; };

inline
std::ostream& operator<< ( std::ostream& os, const tdbscan::Position1d& p1d) {
  return os << std::format("Position1d({})", p1d.value_);
};
} //namespace tdbscan

template <>
struct std::formatter<tdbscan::Position1d> {
  constexpr auto parse(std::format_parse_context & ctx) {return ctx.begin();}

  auto format(const tdbscan::Position1d& p1d, std::format_context &ctx) const {
    return std::format_to(ctx.out(), "[x:{}]", static_cast<double>(p1d));
  };
};


namespace tdbscan {


/// get the distance with a partner object
Position3d::Distance_t
Position3d::distance(const Position3d& rhs) const
		{ return sqrt(pow(xord - rhs.xord, 2) + pow(yord - rhs.yord, 2) + pow(zord - rhs.zord, 2));};

Position3d::Distance_t
Position3d::magnitude() const
		{return sqrt(pow(xord, 2) + pow(yord, 2) + pow(zord,2));};

Position3d::Distance_t Position3d::abs() const
		{return this->magnitude();};

bool Position3d::operator<(const Position3d& rhs) const
		{return magnitude() < rhs.magnitude() || xord < rhs.xord || yord < rhs.yord || zord < rhs.zord ;};

bool Position3d::operator==(const Position3d& rhs) const
		{return xord == rhs.xord && yord == rhs.yord && zord == rhs.zord;};


inline std::ostream& operator<<(std::ostream& os, const Position3d & p3d) {
  return os << std::format("Position3d(x:{}, y:{}, z:{})", p3d.xord, p3d.yord, p3d.zord);
};

}; //namespace tdbscan


template <>
struct std::formatter<tdbscan::Position3d> {
  constexpr auto parse(std::format_parse_context & ctx) {return ctx.begin();}

  auto format(const tdbscan::Position3d& p3d, std::format_context &ctx) const {
    return std::format_to(ctx.out(), "[x:{}, y:{}, z:{}]", p3d.xord, p3d.yord, p3d.zord);
  };
};


#endif //TDBSCAN__COMMON_DEFS_HH
