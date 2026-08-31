#include "tdbscan_algo/tdbscan_algo.h"

#include <iostream>
#include <cstdlib>

#include "tdbscan_algo/common_defs.h"

using namespace std;
using namespace tdbscan;

// define time as a pure positive unit


// define some Limiters
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

// std::set<Blib4d, Blib4d::TimeOrder> construct_blibs() {
std::set<Blib4d> construct_blibs() {
	int many_blibs = 1000;

	std::set<Blib4d> blibs;
	for (int i = 0; i < many_blibs; i++) {
		blibs.insert(Blib4d(rand_pos(), rand_time()));
	}

	return blibs;
}



int main(int argc, char **argv) {
	auto my_algo = construct_algo();

	auto blibs = construct_blibs();

	auto result = my_algo.Process(blibs);

	return 0;
}
