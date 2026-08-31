//
// Created by netsu on 30/08/2026.
//

#ifndef TDBSCAN_COMMON_DEFS_H
#define TDBSCAN_COMMON_DEFS_H

#include <cmath>
#include "tdbscan_algo/base_defs.h"
#include "tdbscan_algo/absblib.h"

// ========================= TIME =======================


/** make a typedef for what is the notion of Time;
* On first principles time is a continuous monotonic increasing variable
*/
class MyTime_t final : tdbscan::Time_t  {
public:
	typedef double TimeDiff_t;
private:
	double value_{0.};
public:
	/// default constructor: for convenience
	MyTime_t() : value_(0.) {};
	MyTime_t(const double value) : value_(value) {};

	inline
	bool
	operator<(const MyTime_t &rhs) const
			{return value_ < rhs.value_;};

	inline
	bool
	operator==(const MyTime_t &rhs) const
			{return value_ == rhs.value_;};

	inline
	TimeDiff_t
	operator-(const MyTime_t &rhs) const
			{return value_ - rhs.value_;};

	//implicit conversion operator for shorthand
	operator double() const { return value_; }
	//assignment operator
	MyTime_t& operator=(const double rhs)
			{value_=rhs; return *this; };

	static double min() {return std::numeric_limits<double>::min();};
	static double max() {return std::numeric_limits<double>::max();};
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

	[[nodiscard]]
	Distance_t
	distance(const Position1d& rhs) const
	{ return value_- rhs.value_;};

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
	operator double() const { return value_; }
	/// assignment operator
	Position1d& operator=(const double rhs)
	{value_=rhs; return *this; };
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
};


// ========================= BLIB =======================
//make a declaration of the Blib
class Blib4d : public tdbscan::AbsBlib<Position3d, MyTime_t> {
public: //type shorthands
	using Ordinate_t = Position3d;
	using Time_t = MyTime_t;

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
	getDistance(const Blib4d& rhs) const
		{return pos.distance(rhs.pos);};

	/// get the time difference
	[[nodiscard]] Time_t::TimeDiff_t
	timeDiff(const Blib4d& other) const
		{return time - other.time;};
public: //comparators
	/// define the lesser-operator
	[[nodiscard]] bool
	operator<(const Blib4d& other) const
		{ return time < other.time || pos < other.pos; };

	[[nodiscard]] bool
	operator==(const Blib4d& other) const
		{ return time == other.time && pos == other.pos; };

	///constructor
	Blib4d(const Position3d pos, const Time_t time) : pos(pos), time(time) {};
};



#endif //TDBSCAN_COMMON_DEFS_H
