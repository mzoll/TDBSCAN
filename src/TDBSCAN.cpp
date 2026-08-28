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
	double value_{0.};

	bool operator<(const MyTime_t &rhs) const {return value_ < rhs.value_;};

	inline
	Timediff_t operator-(const MyTime_t &rhs) const {return value_ - rhs.value_;};
};

// make a definition of a Point in 3d space
class Position3d final : public Position_t {
public:
	double xord, yord, zord;
public:
	/// constructor
	Position3d(const double x, const double y, const double z) : xord(x), yord(y), zord(z) {};

	/// get the distance with a partner object
	[[nodiscard]]
	Distance_t distance(const Position3d& rhs) const { return sqrt(pow(xord - rhs.xord, 2) + pow(yord - rhs.yord, 2) + pow(zord - rhs.zord, 2));};

	[[nodiscard]]
	Distance_t magnitude() const {return sqrt(pow(xord, 2) + pow(yord, 2) + pow(zord,2));};

	[[nodiscard]]
	inline
	Distance_t abs() const {return  this->magnitude();};


	[[nodiscard]]
	bool operator<(const Position3d& rhs) const {return magnitude() < rhs.magnitude() || xord < rhs.xord || yord < rhs.yord || zord < rhs.zord ;};

	[[nodiscard]]
	bool operator==(const Position3d& rhs) const {return xord == rhs.xord && yord == rhs.yord && zord == rhs.zord;};
};


	//make a declaration of the Blib
class Blib4d : public AbsBlib<Position3d> {
private:
	Position3d pos;
	Time_t time;
public:
	[[nodiscard]] virtual Position3d GetPosition() const {return pos;};

	[[nodiscard]] virtual Time_t GetTime() const {return time;};
public: //comparators
	[[nodiscard]] virtual Distance_t GetDistance(const Blib4d& rhs) const {return pos.distance(rhs.pos);};

	/// get the time difference to rhs
	[[nodiscard]] virtual Timediff_t TimeDiff(const Blib4d& other) const {return time - other.time;};

	/// define the lesser-operator
	[[nodiscard]] virtual bool operator<(const Blib4d& other) const {
		return time < other.time || pos < other.pos;
	};

	[[nodiscard]] virtual bool operator==(const Blib4d& other) const {
		return time == other.time && pos == other.pos;
	};

	///constructor
	Blib4d(const Position3d pos, const Time_t time) : pos(pos), time(time) {};
};

// having implemented all this stuff, lets whip up a connector, which is just an distance connector
class DistanceLimiter final : public ConnectorSingle<Blib4d> {
public:
	Distance_t maxDist_;
	DistanceLimiter(const Distance_t maxDistance) : ConnectorSingle("DistConnector"), maxDist_(maxDistance) {};

	bool eval(const Blib4d& lhs, const Blib4d& rhs) const {return lhs.GetDistance(rhs) <= maxDist_;};
};

// make one connector which just connects to max time-diff
class TimeLimiter final : public ConnectorSingle<Blib4d> {
public:
	Timediff_t maxTimediff_;
	explicit TimeLimiter(const Timediff_t maxTimeDiff) : ConnectorSingle("DistConnector"), maxTimediff_(maxTimeDiff) {};

	bool eval(const Blib4d& lhs, const Blib4d& rhs) const {return rhs.TimeDiff(lhs) <= maxTimediff_;};
};

// combine the Connectors into a ConnectorBlock
class LimitingConnector final : public ConnectorBlock<Blib4d> {
};


void construct_algo() {

	auto distLimiter_ = new DistanceLimiter(20.);
	auto timeLimiter_ = new TimeLimiter(20.);
	auto limcon = new LimitingConnector();
	limcon->addConnector(distLimiter_);
	limcon->addConnector(timeLimiter_);


	TDBScan_ParameterSet params;

	params.multiplicity=4;
	params.multiplicityTimeWindow=20;

	TDBScan_Algo<Blib4d> my_algo(params, limcon);
}


double rand_ord() {return rand() * 100-50.;};
Position3d rand_pos() {return Position3d(rand_ord(), rand_ord(), rand_ord())};
Time_t rand_time() {return Time_t(rand() % 10000);};

void construct_blibs() {


	std::set<Blib4d> blibs;
	for (int i = 0; i < 10; i++) {
		blibs.insert(Blib4d(rand_pos(), rand_time()));
	}


}





static void print_hello_world() {
	cout << "Hello world";
}


int main(int argc, char **argv) {
	print_hello_world();
	return 0;
}
