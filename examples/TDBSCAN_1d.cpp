#include "tdbscan_algo/tdbscan_algo.h"

#include <cstdlib>
#include <format>
#include <exception>
#include <fmt/format.h>

#include "tdbscan_algo/common_defs.h"

#include <random>
#include <sstream>

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

/**
 * A twist to the ScalarBlib class that allows to keep track of its source, either Noise or Signal
 */
class SBlibWithTrace final : public ScalarBlib {
public:
  /// denotes the origin of the Blib; either Noise or Signal
	enum Origin {
		UNKNOWN = 0,
		NOISE = 20,
		SIGNAL =100,
	} origin_{UNKNOWN};

  inline
  static std::string origin_tostr(const Origin o) {
    switch (o) {
      case UNKNOWN: return "UNKNOWN";
      case NOISE: return "NOISE";
      case SIGNAL: return "SIGNAL";
      default: throw std::invalid_argument("Value_error");
    }

  }

  ///mark this Blib as to stem either from a Noise or Signal source
	SBlibWithTrace& mark(const Origin o) { origin_ = o; return *this; }

  /// constructor
	SBlibWithTrace( const SBlibWithTrace::Ordinate_t& ord , const SBlibWithTrace::Time_t& t, const Origin o ) : ScalarBlib(ord, t) ,origin_(o) {};

public:
  friend
  std::ostream& operator<< ( std::ostream& outs, const SBlibWithTrace & st);
};

inline
std::ostream& operator<< ( std::ostream& os, const SBlibWithTrace & sblib_trace) {
  return os <<
    SBlibWithTrace::origin_tostr(sblib_trace.origin_) << "::" <<
      static_cast<ScalarBlib>(sblib_trace);
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
	params.emergenceTimeWindow=2;
	params.earlyMergeOverlapRatio= 1.;
	params.lateMergeOverlapRatio= 1.;

	return TDBScan_Algo<SBlibWithTrace>(params, limcon);
}



std::set<SBlibWithTrace>
generate_noise(const double noise_freq, const double width_fields, const double time_duration) {
	std::set<SBlibWithTrace> blibs;
	for (int time_step = 0; time_step < time_duration; time_step++) {
		for (int count_noise = 0; count_noise < noise_freq * width_fields; count_noise++) {
			const double pos = rand_double() * width_fields;
			const double t = rand_double() + time_step;
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

	const auto _box_blibs = generate_moving_box( 5, 2, 0, 1, 50 );
	log_info(std::format("Generated {} BOX blibs", _box_blibs.size()));

	const auto _noise_blibs = generate_noise(0.1, 100, 50.);
  log_info(std::format("Generated {} NOISE blibs", _noise_blibs.size()));
  blibs.insert(_box_blibs.cbegin(), _box_blibs.cend());
	//blibs.insert(_noise_blibs.cbegin(), _noise_blibs.cend());
	return blibs;
}

int main(int argc, char **argv) {
	auto my_algo = construct_algo();

	log_info(std::format("Generate blibs"));
	const auto blibs = gernerate_blibs(100, 50  );

	//take first 3
	std::set<SBlibWithTrace> _blibs;
	auto iter = blibs.begin();
	for (int i = 0; i < 100; i++) {
		_blibs.insert(*iter);
		log_trace( std::ostringstream() << "===sorting==PROBE : {}" << *iter);
		++iter;
	}
  log_info(std::format("Processing nBlibs: {}", _blibs.size()));
	const auto result = my_algo.Process(_blibs);

	log_info(std::format("Generated nClusters: {}", result.size()));

	// for (const auto& c : result) {
	// 	log_info(std::format("Size: {}", c.size()));
	// }
}
