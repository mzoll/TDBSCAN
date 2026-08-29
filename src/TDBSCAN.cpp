#include <iostream>
#include <math.h>
#include <cstdlib>

#include "tdbscan_algo/tdbscan_algo.h"

using namespace std;
using namespace tdbscan;

// define time as a pure positive unit



/** make a typedef for what is the notion of Time;
* On first principles time is a continuous monotonic increasing variable
*/
class MyTime_t final {
public:
	typedef double TimeDiff_t;
private:
	double value_{0.};
public:

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

// make a definition of a Point in 3d space
class Position3d final : public Ordinate_t {
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


	//make a declaration of the Blib
class Blib4d : public AbsBlib<Position3d, MyTime_t> {
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

// having implemented all this stuff, lets whip up a connector, which is just an distance connector
class DistanceLimiter final : public ConnectorSingle<Blib4d> {
public:
	Blib4d::Ordinate_t::Distance_t maxDist_;
	DistanceLimiter(const Blib4d::Ordinate_t::Distance_t maxDistance) : ConnectorSingle("DistConnector"), maxDist_(maxDistance) {};

	bool eval(const Blib4d& lhs, const Blib4d& rhs) const {return lhs.getDistance(rhs) <= maxDist_;};
};

// make one connector which just connects to max time-diff
class TimeLimiter final : public ConnectorSingle<Blib4d> {
public:
	Blib4d::Time_t::TimeDiff_t maxTimediff_;
	explicit TimeLimiter(const Blib4d::Time_t::TimeDiff_t maxTimeDiff) : ConnectorSingle("DistConnector"), maxTimediff_(maxTimeDiff) {};

	bool eval(const Blib4d& lhs, const Blib4d& rhs) const {return rhs.timeDiff(lhs) <= maxTimediff_;};
};

// combine the Connectors into a ConnectorBlock
class LimitingConnector final : public ConnectorBlock<Blib4d> {
};


TDBScan_Algo<Blib4d> construct_algo() {

	auto distLimiter_ = new DistanceLimiter(20.);
	auto timeLimiter_ = new TimeLimiter(20.);
	auto limcon = new LimitingConnector();
	limcon->addConnector(distLimiter_);
	limcon->addConnector(timeLimiter_);


	TDBScan_Algo<Blib4d>::TDBScan_ParameterSet params;

	params.multiplicity=4;
	params.multiplicityTimeWindow=20;

	return TDBScan_Algo<Blib4d>(params, limcon);
}


double rand_ord() {return rand() * 100-50.;};
Position3d rand_pos() {return Position3d(rand_ord(), rand_ord(), rand_ord());};
MyTime_t rand_time() {return MyTime_t(rand() % 10000);};

std::set<Blib4d> construct_blibs() {
	int many_blibs = 1000;

	std::set<Blib4d> blibs;
	for (int i = 0; i < many_blibs; i++) {
		blibs.insert(Blib4d(rand_pos(), rand_time()));
	}

	return blibs;
}


static void print_hello_world() {
	cout << "Hello world";
}


int main(int argc, char **argv) {
	auto my_algo = construct_algo();

	auto blibs = construct_blibs();

	auto result = my_algo.Process(blibs);

	return 0;
}
