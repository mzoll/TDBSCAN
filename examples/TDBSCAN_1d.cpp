#include "tdbscan_algo/tdbscan_algo.h"

#include <cstdlib>
#include <format>

#include "tdbscan_algo/common_defs.h"

#include <random>

using namespace std;
using namespace tdbscan;

// define time as a pure positive unit

double rand_double() {
	double lower_bound = 0.;
	double upper_bound = 1.;
	static std::uniform_real_distribution<double> unif(lower_bound,upper_bound);
	static std::default_random_engine re;
	return unif(re);
}



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

	auto distLimiter_ = new DistanceLimiter(4.);
	auto timeLimiter_ = new TimeLimiter(2.);
	auto limcon = new LimitingConnector();
	limcon->addConnector(distLimiter_);
	limcon->addConnector(timeLimiter_);

	TDBScan_Algo<ScalarBlib>::TDBScan_ParameterSet params;

	params.multiplicity=4;
	params.multiplicityTimeWindow=2;
	params.earlyMergeOverlapRatio= 1.;
	params.lateMergeOverlapRatio= 1.;

	return TDBScan_Algo<ScalarBlib>(params, limcon);
}


std::set<ScalarBlib>
generate_noise(const double noise_freq, const double width_fields, const double time_duration) {
	std::set<ScalarBlib> blibs;
	for (int time_step = 0; time_step < time_duration * noise_freq; time_step++) {
		const auto pos = rand_double() * width_fields;
		const auto t = rand_double() * time_duration;
		log_trace(std::format("ONE: {}", time_step));
		blibs.insert(ScalarBlib({pos}, t));
	}

	return blibs;
}

std::set<ScalarBlib>
generate_moving_box(const double box_size, const double inerta, const double start_pos, const double brightness, const double time_duration) {
	std::set<ScalarBlib> blibs;
	const auto _brightness_cal = brightness * box_size;

	for (int time_step = 0; time_step < time_duration; time_step++) {
		for (int j = 0; j < _brightness_cal; j++) {
			const auto t = rand_double() + time_step;
			const auto box_ledge_pos = time_step * inerta + start_pos - box_size /2.;
			const auto pos = rand_double() * box_size + box_ledge_pos;
			blibs.insert(ScalarBlib({pos}, t));
		}
	}
	return blibs;
}

std::set<ScalarBlib>
gernerate_blibs( const double width_fields, const double time_duration ) {
	std::set<ScalarBlib> blibs;

	log_info(std::format("Generate BOX blibs"));
	const auto _box_blibs = generate_moving_box( 5, 2, 0, 1, 50 );
	log_info(std::format("Generate NOISE blibs"));
	const auto _noise_blibs = generate_noise(0.1, 100, 50.);
	blibs.insert(_box_blibs.cbegin(), _box_blibs.cend());
	blibs.insert(_noise_blibs.cbegin(), _noise_blibs.cend());
	return blibs;
}

int main(int argc, char **argv) {
	auto my_algo = construct_algo();

	log_info(std::format("Generate blibs"));
	const auto blibs = gernerate_blibs(100, 50  );
	log_info(std::format("Processing nBlibs: {}", blibs.size()));
	//take first 3
	std::set<ScalarBlib> _blibs;
	auto iter = blibs.begin();
	for (int i = 0; i < 3; i++) {
		_blibs.insert(*iter);
		++iter;
	}
	const auto result = my_algo.Process(_blibs);

	log_info(std::format("Generated nClusters: {}", result.size()));

	// for (const auto& c : result) {
	// 	log_info(std::format("Size: {}", c.size()));
	// }
}
