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



class SBlibWithTrace final : public ScalarBlib {
public:
	enum Origin {
		UNKNOWN = 0,
		NOISE = 20,
		SIGNAL =100,
	} origin_{UNKNOWN};

	SBlibWithTrace& mark(const Origin o) { origin_ = o; return *this; }

	SBlibWithTrace( const SBlibWithTrace::Ordinate_t& ord , const SBlibWithTrace::Time_t& t, const Origin o ) : ScalarBlib(ord, t) ,origin_(o) {};
};



// define some Limiters
class DistanceLimiter final : public ConnectorSingle<SBlibWithTrace> {
public:
	SBlibWithTrace::Ordinate_t::Distance_t maxDist_;
	DistanceLimiter(const SBlibWithTrace::Ordinate_t::Distance_t maxDistance) : ConnectorSingle("DistConnector"), maxDist_(maxDistance) {};

	bool eval(const SBlibWithTrace& lhs, const SBlibWithTrace& rhs) const {return lhs.getDistance(rhs) <= maxDist_;};
};

// make one connector which just connects to max time-diff
class TimeLimiter final : public ConnectorSingle<SBlibWithTrace> {
public:
	SBlibWithTrace::Time_t::TimeDiff_t maxTimediff_;
	explicit TimeLimiter(const SBlibWithTrace::Time_t::TimeDiff_t maxTimeDiff) : ConnectorSingle("DistConnector"), maxTimediff_(maxTimeDiff) {};

	bool eval(const SBlibWithTrace& lhs, const SBlibWithTrace& rhs) const {return rhs.timeDiff(lhs) <= maxTimediff_;};
};

// combine the Connectors into a ConnectorBlock
class LimitingConnector final : public ConnectorBlock<SBlibWithTrace> {
};





TDBScan_Algo<SBlibWithTrace> construct_algo() {

	auto distLimiter_ = new DistanceLimiter(4.);
	auto timeLimiter_ = new TimeLimiter(2.);
	auto limcon = new LimitingConnector();
	limcon->addConnector(distLimiter_);
	limcon->addConnector(timeLimiter_);

	TDBScan_Algo<SBlibWithTrace>::TDBScan_ParameterSet params;

	params.multiplicity=4;
	params.multiplicityTimeWindow=2;
	params.earlyMergeOverlapRatio= 1.;
	params.lateMergeOverlapRatio= 1.;

	return TDBScan_Algo<SBlibWithTrace>(params, limcon);
}



std::set<SBlibWithTrace>
generate_noise(const double noise_freq, const double width_fields, const double time_duration) {
	std::set<SBlibWithTrace> blibs;
	for (int time_step = 0; time_step < time_duration; time_step++) {
		for (int count_noise = 0; count_noise < noise_freq * width_fields; count_noise++) {
			const auto pos = rand_double() * width_fields;
			const auto t = rand_double() + time_step;
			blibs.insert(SBlibWithTrace({pos}, t, SBlibWithTrace::NOISE));
		}
	}
	return blibs;
}


std::set<SBlibWithTrace>
generate_moving_box(const double box_size, const double inerta, const double start_pos, const double brightness, const double time_duration) {
	std::set<SBlibWithTrace> blibs;
	const auto _brightness_cal = brightness * box_size;

	for (int time_step = 0; time_step < time_duration; time_step++) {
		for (int j = 0; j < _brightness_cal; j++) {
			const auto t = rand_double() + time_step;
			const auto box_ledge_pos = time_step * inerta + start_pos - box_size /2.;
			const auto pos = rand_double() * box_size + box_ledge_pos;
			blibs.insert(SBlibWithTrace({pos}, t, SBlibWithTrace::SIGNAL));
		}
	}
	return blibs;
}

std::set<SBlibWithTrace>
gernerate_blibs( const double width_fields, const double time_duration ) {
	std::set<SBlibWithTrace> blibs;

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
	std::set<SBlibWithTrace> _blibs;
	auto iter = blibs.begin();
	for (int i = 0; i < 10; i++) {
		_blibs.insert(*iter);
		log_trace(std::format("===sorting==PROBE : o:{} t:{}", double(iter->getOrdinate()), double(iter->getTime() )));
		++iter;
	}
	const auto result = my_algo.Process(_blibs);

	log_info(std::format("Generated nClusters: {}", result.size()));

	// for (const auto& c : result) {
	// 	log_info(std::format("Size: {}", c.size()));
	// }
}
