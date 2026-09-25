//
// Created by mzoll on 01/08/2026.
//

#include "tdbscan_algo/tdbscan_algo.h"

#include <cstdlib>
#include <format>
#include <exception>
#include <cmath>
#include <fmt/format.h>

#include "tdbscan_algo/common_defs.h"

#include <random>
#include <sstream>

using namespace std;
using namespace tdbscan;


// create a random double on
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

  ///convenience to convert `Origin` into a string for human reading
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

  /**
   * Fully qualified constructor
   * @param ord the Ordinate
   * @param t time of the blib
   * @param o the origin of the Blib, aka SIGNAL, or NOISE
   */
  SBlibWithTrace(
	  const SBlibWithTrace::Ordinate_t& ord,
	  const SBlibWithTrace::Time_t& t,
	  const Origin o ) :
  ScalarBlib(ord, t) ,origin_(o) {};

public:
  friend
  std::ostream& operator<< ( std::ostream& outs, const SBlibWithTrace & st);
};

/// custom ostream output; nicely formats the class as a string
inline
std::ostream& operator<< ( std::ostream& os, const SBlibWithTrace & sblib_trace) {
  return os <<
    SBlibWithTrace::origin_tostr(sblib_trace.origin_) << "::" <<
      static_cast<ScalarBlib>(sblib_trace);
};

/// custom formatter for SBlibWithTrace. Used in conjunction with `std::format` or `fmt::format`
template <>
struct std::formatter<SBlibWithTrace> {
  constexpr auto parse(std::format_parse_context & ctx) {return ctx.begin();}

  auto format(const SBlibWithTrace& sb, std::format_context &ctx) const {
    return std::format_to(ctx.out(), "[ord:{}, time:{}]::{}",
      static_cast<double>(sb.getOrdinate()),
      static_cast<double>(sb.getTime()),
      SBlibWithTrace::origin_tostr(sb.origin_));
  };
};


// ===================== Define some Limiters ====================
// We create two simple limiters: One to limit distance and one to Limit time between Blibs
// Both of them are Connectors and compare one Blib to another.
// The only thing we need to implement is the the `bool eval(Blib rhs, Blib rhs)´ function.
// In the end we combine both Limiters into a single `ConnectorBlock`, which than can be used in the TDBscan-algorithm
// ========================================================
class DistanceLimiter final : public ConnectorSingle<SBlibWithTrace> {
public:
	SBlibWithTrace::Ordinate_t::Distance_t maxDist_;
	DistanceLimiter(const SBlibWithTrace::Ordinate_t::Distance_t maxDistance) :
  ConnectorSingle("DistConnector"), maxDist_(maxDistance) {};

	bool eval(const SBlibWithTrace& lhs, const SBlibWithTrace& rhs) const {
	  return abs(lhs.distanceTo(rhs)) <= maxDist_;
	};
};

// make one connector which just connects to max time-diff; this is time ordered and thereby is positive one-sided
class TimeLimiter final : public ConnectorSingle<SBlibWithTrace> {
public:
	SBlibWithTrace::Time_t::TimeDiff_t maxTimediff_;
	explicit TimeLimiter(const SBlibWithTrace::Time_t::TimeDiff_t maxTimeDiff) :
  ConnectorSingle("DistConnector"), maxTimediff_(maxTimeDiff) {};

	bool eval(const SBlibWithTrace& lhs, const SBlibWithTrace& rhs) const {
	  return rhs.timeTo(lhs) <= maxTimediff_;
	};
};

// combine the Connectors into a ConnectorBlock
class LimitingConnector final : public ConnectorBlock<SBlibWithTrace> {
};



/* ============================ Constructing the TDBscan algorithm instance ==================
 * the main algorithm is constructed by creating instances of the limiters defined in the previous set,
 * and providing a set of Algorithm parameters
 */
TDBScan_Algo<SBlibWithTrace> construct_algo(const double distance_lim, const double time_lim) {

	auto distLimiter_ = new DistanceLimiter(distance_lim);
	auto timeLimiter_ = new TimeLimiter(time_lim);
	auto limcon = new LimitingConnector();
	limcon->addConnector(distLimiter_);
	limcon->addConnector(timeLimiter_);

	TDBScan_Algo<SBlibWithTrace>::TDBScan_ParameterSet params;

	params.multiplicity=4;
	params.multiplicityTimeWindow=2.; //internally converted into SBlibWithTrace<..., tTime>::tTime::Time_t
	params.emergenceTimeWindow=2.;//internally converted into SBlibWithTrace<..., tTime>::tTime::Time_t
	params.earlyMergeOverlapRatio= 1.; //ZeroOne: its a ratio
	params.lateMergeOverlapRatio= 1.; //ZeroOne: its a ratio

	return TDBScan_Algo<SBlibWithTrace>(params, limcon);
}


/* ========================= Create the scenario =======================
 * For a demonstration and testing scenario Blibs need to be generated which stem from two distinguishable sources:
 * Signal and Noise.
 * In Our case of a 1d- grid, the signal will be a represented by a box of finite size moving sideways
 * with a given inertia generating blibs according to a brightness. Noise on the other hand will be random noise
 * uniformly generated on the fields of that grid.
 */
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


/**
 * generate blibs from a box moving over a 1d space left to right
 *
 * example: box_size: 4, inertia: 2, ledge_start_pos: 0, brightness: 0.5, time_duration 5
 * --[xx  ]-------------------
 * --.--[ x x]----------------
 * --.----[x  x]--------------
 * --.------[ xx ]------------
 * --.--------[  xx]----------
 *
 * The box will be reflected on the right and left side
 *
 * @param box_size
 * @param inerta
 * @param start_pos,
 * @param brightness a measure of the signal frequency in one unit-volume of the box
 * @param field_size,
 * @param time_duration
 * @return
 */
std::set<SBlibWithTrace>
generate_moving_box(
  const double box_size,
  const double inertia,
  const double start_pos,
  const double brightness,
  const double field_size,
  const double time_duration) {
	std::set<SBlibWithTrace> blibs;

  const int n_blibs = static_cast<int>(time_duration * box_size * brightness);

  for (int i_blib = 0; i_blib < n_blibs; i_blib++) {
    const auto t = rand_double() * time_duration;
    const auto center_pos_ind = (t*inertia + start_pos) / field_size;

    const double center_pos = field_size * (int(center_pos_ind) % 2 == 0 ?  center_pos_ind - std::floor(center_pos_ind) : 1. - (center_pos_ind - std::floor(center_pos_ind)));

    const double blib_pos_inbox = rand_double() * box_size - box_size/2.;

    const auto blib_pos = center_pos + blib_pos_inbox;

    if (blib_pos<0. || blib_pos > field_size)
      continue;

    blibs.insert(SBlibWithTrace({blib_pos}, t, SBlibWithTrace::SIGNAL));
  }

	return blibs;
}


/**
 * Generate Blibs for our scenario
 *
 * A bright box moves left to right in
 * @param time_duration
 * @param width_fields
 * @param brightness
 * @param noise_contamination
 * @return
 */
std::set<SBlibWithTrace>
gernerate_blibs( const double time_duration=50, const double width_fields = 100, const double brightness= 1., const double noise_contamination = 0.1) {
	std::set<SBlibWithTrace> blibs;

	const auto _box_blibs = generate_moving_box(5., 1., 0., brightness, width_fields, time_duration);
	LOG_INFO("Generated {} BOX blibs", _box_blibs.size());
  blibs.insert(_box_blibs.cbegin(), _box_blibs.cend());

	const auto _noise_blibs = generate_noise(noise_contamination*brightness, width_fields, time_duration);
  LOG_INFO("Generated {} NOISE blibs", _noise_blibs.size());
	//blibs.insert(_noise_blibs.cbegin(), _noise_blibs.cend());
	return blibs;
}


/**
 * Helper function: calculate the purity, Signal over Noise ratio, of this Blib sample
 * @param blibs the blibs to quantify
 * @return the purity as '#signal / #all'.
 */
double calculate_signal_purity(const std::set<SBlibWithTrace>& blibs) {
  int _signal_count = 0;
  int _noise_count = 0;
  for (const auto& b: blibs) {
    switch (b.origin_) {
      case SBlibWithTrace::SIGNAL:
        _signal_count++;
        break;
      case SBlibWithTrace::NOISE:
        _noise_count++;
        break;
      case SBlibWithTrace::UNKNOWN:
        break;
    }
  }
  return static_cast<double>(_signal_count)/blibs.size();
}



int main(int argc, char **argv) {
	auto my_algo = construct_algo(2., 0.5);

	LOG_INFO("Generate blibs");
	const auto blibs = gernerate_blibs(50 ,100, 3, 0.1 );

  LOG_INFO("Processing nBlibs: {} (purity {:.3f})", blibs.size(), calculate_signal_purity(blibs));
  //take first 3
  std::set<SBlibWithTrace> _blibs;
  auto iter = blibs.begin();
  for (int i = 0; i < 20; i++) {
    _blibs.insert(*iter);
    LOG_TRACE( "Sample : {}", *iter);
    ++iter;
  }


	const auto result = my_algo.Process(_blibs);

	LOG_INFO("Generated nClusters: {}", result.size());

  auto r_citer = result.cbegin();
  for (int i =0; i< 5; i++) {
    if (r_citer == result.cend())
      break;
    LOG_INFO("Cluster {} size: {} (purity {:.3f})", i, r_citer->size(), calculate_signal_purity(*r_citer));
    for (const auto& b : *r_citer) {
      LOG_INFO("das {}", b)
    }

    r_citer++;
  }

	// for (const auto& c : result) {
	// 	LOG_INFO(std::format("Size: {}", c.size()));
	// }
}
