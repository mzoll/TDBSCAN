//
// Created by mzoll on 30/08/2026.
//

#ifndef TDBSCAN_COMMON_DEFS_H
#define TDBSCAN_COMMON_DEFS_H

#include <cmath>
#include <ostream>
#include <format>
#include "tdbscan_algo/base_defs.h"
#include "tdbscan_algo/absblib.h"

// ========================= TIME =======================

/** make a typedef for what is the notion of Time;
* On first principles time is a continuous monotonic increasing variable
*/
class ScalarTime_t final : tdbscan::Time_t  {
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
	operator<(const ScalarTime_t &rhs) const
			{return value_ < rhs.value_;};

	inline
	bool
	operator==(const ScalarTime_t &rhs) const
			{return value_ == rhs.value_;};

	inline
	TimeDiff_t
	operator-(const ScalarTime_t &rhs) const
			{return value_ - rhs.value_;};

	//implicit conversion operator for shorthand
	operator double() const { return value_; }
	//assignment operator
	ScalarTime_t& operator=(const double rhs)
			{value_=rhs; return *this; };

	static constexpr double min() {return -std::numeric_limits<double>::infinity();};
	static constexpr double max() {return std::numeric_limits<double>::infinity();};

private:
	friend
	std::ostream& operator<< ( std::ostream& os, const ScalarTime_t & st);
};

inline
std::ostream& operator<< ( std::ostream& os, const ScalarTime_t & st) {
	return os << std::format("ScalarTime({})", st.value_);
};


template <>
struct std::formatter<ScalarTime_t> {
  constexpr auto parse(std::format_parse_context & ctx) {return ctx.begin();}

  auto format(const ScalarTime_t& st, std::format_context &ctx) const {
    return std::format_to(ctx.out(), "t_{}", static_cast<double>(st));
  };
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
	distance(const Position1d& rhs) const
	{ return rhs.value_ - value_;};

	[[nodiscard]]
	Distance_t magnitude() const
	{return value_;};

	[[nodiscard]]
	inline
	Distance_t abs() const
	{return this->magnitude();};

	bool
	operator<(const Position1d &rhs) const
	{return value_ < rhs.value_;};

	bool
	operator==(const Position1d &rhs) const
	{return value_ == rhs.value_;};
public: //convenience
	///implicit conversion operator for shorthand
  explicit operator double() const { return value_; }
	/// assignment operator
	Position1d& operator=(const double rhs)
	  {value_=rhs; return *this; };

private:
  friend
  std::ostream& operator<< ( std::ostream& os, const Position1d & p1d);
};


inline std::ostream& operator<< ( std::ostream& os, const Position1d & p1d) {
  return os << std::format("Position1d({})", p1d.value_);
};

template <>
struct std::formatter<Position1d> {
  constexpr auto parse(std::format_parse_context & ctx) {return ctx.begin();}

  auto format(const Position1d& p1d, std::format_context &ctx) const {
    return std::format_to(ctx.out(), "[x:{}]", static_cast<double>(p1d));
  };
};



// make a definition of a Point in 3d space
class Position3d final : public tdbscan::Ordinate_t {
public:
	typedef double Distance_t;
public:
	double xord, yord, zord;
public:
	/// constructor
	Position3d(const double x, const double y, const double z) :
			xord(x), yord(y), zord(z) {};

	/// get the distance with a partner object
	[[nodiscard]]
	Distance_t
	distance(const Position3d& rhs) const
			{ return sqrt(pow(xord - rhs.xord, 2) + pow(yord - rhs.yord, 2) + pow(zord - rhs.zord, 2));};

	[[nodiscard]]
	Distance_t magnitude() const
			{return sqrt(pow(xord, 2) + pow(yord, 2) + pow(zord,2));};

	[[nodiscard]]
	inline
	Distance_t abs() const
			{return this->magnitude();};

	[[nodiscard]]
	bool operator<(const Position3d& rhs) const
			{return magnitude() < rhs.magnitude() || xord < rhs.xord || yord < rhs.yord || zord < rhs.zord ;};

	[[nodiscard]]
	bool operator==(const Position3d& rhs) const
			{return xord == rhs.xord && yord == rhs.yord && zord == rhs.zord;};

private:
  friend
  std::ostream& operator<< ( std::ostream& os, const Position3d & p3d);
};

inline std::ostream& operator<< ( std::ostream& os, const Position3d & p3d) {
  return os << std::format("Position3d(x:{}, y:{}, z:{})", p3d.xord, p3d.yord, p3d.zord);
};

template <>
struct std::formatter<Position3d> {
  constexpr auto parse(std::format_parse_context & ctx) {return ctx.begin();}

  auto format(const Position3d& p3d, std::format_context &ctx) const {
    return std::format_to(ctx.out(), "[x:{}, y:{}, z:{}]", p3d.xord, p3d.yord, p3d.zord);
  };
};



// ========================= BLIB =======================
//make a declaration of the Blib
class Blib3d : public tdbscan::AbsBlib<Position3d, ScalarTime_t> {
public: //type shorthands
	using Ordinate_t = Position3d;
	using Time_t = ScalarTime_t;

private:
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
	getDistance(const Blib3d& rhs) const
		{return pos.distance(rhs.pos);};

	/// get the time difference
	[[nodiscard]] Time_t::TimeDiff_t
	timeDiff(const Blib3d& other) const
		{return time - other.time;};
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

inline std::ostream& operator<< ( std::ostream& os, const Blib3d & b4d) {
  return os << std::format("Blib3d(ord:{}, time:{})", b4d.getOrdinate(), b4d.getTime());
};

// template <>
// struct std::formatter<Blib4d> {
//   constexpr auto parse(std::format_parse_context & ctx) {return ctx.begin();}
//
//   auto format(const Blib4d& sb, std::format_context &ctx) const {
//     return std::format_to(ctx.out(), "[ord:{}, time:{}]", sb.getOrdinate()), static_cast<double>(sb.getTime());
//   };
// };

/**
 * A Blib that has a scalar Ordinate
 */
class ScalarBlib : public tdbscan::AbsBlib<Position1d, ScalarTime_t> {
public: //type shorthands
	using Ordinate_t = Position1d;
	using Time_t = ScalarTime_t;

private:
	Ordinate_t pos;
	Time_t time;
public:
	[[nodiscard]] Position1d
	getOrdinate() const
	{return pos;};

	[[nodiscard]] Time_t
	getTime() const
		{return time;};

  ///
	[[nodiscard]] Ordinate_t::Distance_t
	distanceTo(const ScalarBlib& rhs) const
		{return pos.distance(rhs.pos);};

	/// get the time difference
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

inline
std::ostream& operator << ( std::ostream& os, const ScalarBlib & sb) {
	return os << std::format("[ord: {}, time: {}]", double(sb.getOrdinate()), double(sb.getTime()));
};

template <>
struct std::formatter<ScalarBlib> {
  constexpr auto parse(std::format_parse_context & ctx) {return ctx.begin();}

  auto format(const ScalarBlib& sb, std::format_context &ctx) const {
    return std::format_to(ctx.out(), "[ord:{}, time:{}]", static_cast<double>(sb.getOrdinate()), static_cast<double>(sb.getTime()) );
  };
};


#endif //TDBSCAN_COMMON_DEFS_H
