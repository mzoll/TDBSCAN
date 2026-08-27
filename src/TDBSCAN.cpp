#include <iostream>
#include <math.h>

#include "tdbscan_algo/tdbscan_algo.h"

using namespace std;
using namespace tdbscan;


// make a definition of a Point in 3d space

class Position3d : public Position_t {
public:
	double xord, yord, zord;
public:
	/// constructor
	Position3d(const double x, const double y, const double z) : xord(x), yord(y), zord(z) {};

	/// get the distance with a partner object
	[[nodiscard]]
	Distance_t distance(const Position3d& rhs) const { return sqrt(pow(xord - rhs.xord, 2) + pow(yord - rhs.yord, 2) + pow(zord - rhs.zord, 2));};
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
	[[nodiscard]] virtual Timediff_t TimeDiff(const Blib4d& rhs) const {return time - rhs.time;};
};

// having implemented all this stuff, lets whip up a connector, which is just an distance connector
class DistanceLimiter : public ConnectorSingle<Blib4d> {
public:
	Distance_t maxDist_;
	DistanceLimiter(const Distance_t maxDistance) : ConnectorSingle("DistConnector"), maxDist_(maxDistance) {};

	bool eval(const Blib4d& lhs, const Blib4d& rhs) const {return lhs.GetDistance(rhs) <= maxDist_;};
};

// make one connector which just connects to max time-diff
class TimeLimiter : public ConnectorSingle<Blib4d> {
public:
	Timediff_t maxTimediff_;
	TimeLimiter(const Timediff_t maxTimeDiff) : ConnectorSingle("DistConnector"), maxTimediff_(maxTimeDiff) {};

	bool eval(const Blib4d& lhs, const Blib4d& rhs) const {return rhs.TimeDiff(lhs) <= maxTimediff_;};
};

// combine the Connectors into a ConnectorBlock
class LimitingConnector : public ConnectorBlock<Blib4d> {
	DistanceLimiter distLimiter_{20.};
	TimeLimiter timeLimiter_{20.};

	LimitingConnector() {
		addConnector(&distLimiter_);
		addConnector(&timeLimiter_);
	};


};






};






int main(int argc, char **argv) {
	cout << "Hello world";
	return 0;
}
