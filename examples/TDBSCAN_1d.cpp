#include "tdbscan_algo/tdbscan_algo.h"

#include <iostream>
#include <cstdlib>

#include "tdbscan_algo/common_defs.h"

using namespace std;
using namespace tdbscan;

// define time as a pure positive unit


// define some Limiters
class DistanceLimiter final : public ConnectorSingle<ScalarBlib> {
public:
	ScalarBlib::Ordinate_t::Distance_t maxDist_;
	DistanceLimiter(const ScalarBlib::Ordinate_t::Distance_t maxDistance) : ConnectorSingle("DistConnector"), maxDist_(maxDistance) {};

	bool eval(const ScalarBlib& lhs, const ScalarBlib& rhs) const {return lhs.getDistance(rhs) <= maxDist_;};
};

// make one connector which just connects to max time-diff
class TimeLimiter final : public ConnectorSingle<ScalarBlib> {
public:
	ScalarBlib::Time_t::TimeDiff_t maxTimediff_;
	explicit TimeLimiter(const ScalarBlib::Time_t::TimeDiff_t maxTimeDiff) : ConnectorSingle("DistConnector"), maxTimediff_(maxTimeDiff) {};

	bool eval(const ScalarBlib& lhs, const ScalarBlib& rhs) const {return rhs.timeDiff(lhs) <= maxTimediff_;};
};

// combine the Connectors into a ConnectorBlock
class LimitingConnector final : public ConnectorBlock<ScalarBlib> {
};


TDBScan_Algo<ScalarBlib> construct_algo() {

	auto distLimiter_ = new DistanceLimiter(20.);
	auto timeLimiter_ = new TimeLimiter(20.);
	auto limcon = new LimitingConnector();
	limcon->addConnector(distLimiter_);
	limcon->addConnector(timeLimiter_);

	TDBScan_Algo<ScalarBlib>::TDBScan_ParameterSet params;

	params.multiplicity=4;
	params.multiplicityTimeWindow=20;

	return TDBScan_Algo<ScalarBlib>(params, limcon);
}


std::set<ScalarBlib>
generate_noise(const double noise_freq, const double width_fields, const double time_duration) {
	std::set<ScalarBlib> blibs;
	for (int i = 0; i < time_duration * noise_freq; i++) {
		pos = rand() * width_fields;
		t = rand() * time_duration;
		blibs.push_back(ScalarBlib({pos}, t));
	}

	return blibs;
}

std::set<ScalarBlib>
generate_moving_box(const double box_size, const double inerta, const double start_pos, const double brightness, const double time_duration) {
	std::set<ScalarBlib> blibs;
	const auto _brightness_cal = brightness * box_size;

	for (int time_step = 0; time_step < time_duration; i++) {
		for (int j = 0; j < _brightness_cal; j++) {
			const auto t = rand() + time_step;
			const auto box_ledge_pos = time_step * inerta + start_pos - box_size /2.;
			const auto pos = rand() * _box_size + box_ledge_pos;
			blibs.push_back(ScalarBlib({pos}, t));
		}
	}
	return blibs;
}

gernerate_blibs( const double width_fields, const double time_duration ) {
	std::set<ScalarBlib> blibs;

	const auto _box_blibs = generate_moving_box( 5, 2, 0, 1, 50 );
	const auto _noise_blibs = generate_noise(0.1, 100, 50.);

	blibs.insert(_box_blibs.cbegin(), _box_blibs.cend());
	blibs.insert(_noise_blibs.cbegin(), _noise_blibs.cend());
}



int main(int argc, char **argv) {
	auto my_algo = construct_algo();

	auto blibs = construct_blibs();

	auto result = my_algo.Process(blibs);

	return 0;
}
