//
// Created by mzoll on 01/08/2026.
//

#include <functional>

#include <random>
#include <utility>
#include "external/common_clib/interrupt.h"

#include "tdbscan_algo/tdbscan_algo.h"
#include "external/common_clib/stopwatch.h"
#include "tdbscan_algo/common_blibs.h"
#include "tdbscan_algo/tdbscan_async.h"

using namespace std;
using namespace tdbscan;
using namespace common_clib;

using namespace std::chrono_literals;



/**
 * Limits the Distance within which ScalarBlibs can causally connect
 */
class DistanceLimiter : public ConnectorSingle<ScalarBlib> {
protected:
  ScalarBlib::Ordinate_t::Distance_t maxDist_;

public:
  DistanceLimiter(const ScalarBlib::Ordinate_t::Distance_t maxDistance) : ConnectorSingle("DistanceLimiter"),
                                                                               maxDist_(maxDistance) {};

  [[nodiscard]] bool eval(const ScalarBlib &lhs, const ScalarBlib &rhs) const override {
    return fabs(lhs.distanceTo(rhs)) <= maxDist_;
  };
};


/**
 * Limits the Distance within which Blibs causally connect
 */
class TimeLimiter : public ConnectorSingle<ScalarBlib> {
private:
  ScalarBlib::Time_t::TimeDiff_t maxTimediff_;

public:
  TimeLimiter(const ScalarBlib::Time_t::TimeDiff_t maxTimeDiff)
  : ConnectorSingle("TimeLimiter"), maxTimediff_(maxTimeDiff) {};

  [[nodiscard]] bool eval(const ScalarBlib &lhs, const ScalarBlib &rhs) const final {
    return fabs(lhs.timeTo(rhs)) <= maxTimediff_;
  }
};


/// combine the Connectors into a ConnectorBlock
class LimitingConnector final : public ConnectorAssembly_AND<ScalarBlib> {
  const DistanceLimiter *const distance_limiter_;
  const TimeLimiter *const time_limiter_;

public:
  LimitingConnector(const double distance_lim, const double time_lim)
  : distance_limiter_(
    new DistanceLimiter(distance_lim)),
    time_limiter_(new TimeLimiter(time_lim)) {
    addConnector(distance_limiter_);
    addConnector(time_limiter_);
  };

  ~LimitingConnector() override {
    delete distance_limiter_;
    delete time_limiter_;
  };
};

/* ============================ Constructing the TDBscan algorithm instance ==================
 * the main algorithm is constructed by creating instances of the limiters defined in the previous set,
 * and providing a set of Algorithm parameters
 */
TDBScan_AsyncMachine<ScalarBlib> construct_algo(const double distance_lim, const double time_lim) {
  auto limcon = new LimitingConnector(distance_lim, time_lim);

  TDBScan_Algo<ScalarBlib>::TDBScan_ParameterSet params;

  params.multiplicity = 4;
  params.multiplicityTimeWindow = 2.; //internally converted into SBlibWithTrace<..., tTime>::tTime::Time_t
  params.emergenceTimeWindow = 2.; //internally converted into SBlibWithTrace<..., tTime>::tTime::Time_t
  params.earlyMergeMultiplicityRatio = 1.; //ZeroOne: its a ratio
  params.lateMergeOverlapRatio = 1.; //ZeroOne: its a ratio

  return TDBScan_AsyncMachine<ScalarBlib>(params, limcon);
}



class threaded {
  mutable std::thread* thread_;
  bool should_run_{false};
  const int spurious_wakeup_time_ms_;

  mutable struct Telemetry {
    unsigned int n_processed{0};
    unsigned int n_paused{0};
  } telemetry_;


protected:
  virtual void worktask() const = 0;
private:
  void worktask_with_telemetry() const { this->worktask() ; ++telemetry_.n_processed; };

public:
  threaded(const int spurious_wakeup_time_ms = 1000) : spurious_wakeup_time_ms_(spurious_wakeup_time_ms) {};
  ~threaded() {stop();}

  void start() {
    should_run_ = true;
    thread_ = new std::thread(&threaded::worktask_with_telemetry, this);
  };

  void stop() {
    should_run_ = false;
    if (thread_) {
      thread_->join();
      delete thread_;
      thread_=nullptr;
    }
  };

};


class Generator {
public :
  struct generator_exhausted : std::exception {};
private:
  // create a random double on
  static double rand_double() {
    double lower_bound = 0.;
    double upper_bound = 1.;
    static std::uniform_real_distribution<double> unif(lower_bound, upper_bound);
    static std::default_random_engine re;
    return unif(re);
  }

  static std::set<ScalarBlib>
  generate_noise(const double noise_freq, const double width_fields, const double time_duration) {
    std::set<ScalarBlib> blibs;
    for (int time_step = 0; time_step < time_duration; time_step++) {
      for (int count_noise = 0; count_noise < noise_freq * width_fields; count_noise++) {
        const double pos = rand_double() * width_fields;
        const double t = rand_double() + time_step;
        blibs.insert(ScalarBlib({pos}, t));
      }
    }
    return blibs;
  }
private:
  mutable set<ScalarBlib> blib_pool_;

public:
  Generator(const double noise_freq, const double width_fields, const double time_duration)
  : blib_pool_(generate_noise(noise_freq, width_fields, time_duration)) {};
public:
  ScalarBlib next() const {
    while (!blib_pool_.empty()) {
      const auto blib = *blib_pool_.cbegin();
      blib_pool_.erase(blib_pool_.cbegin());
      return blib;
    }
    throw generator_exhausted();
  };
};


//========================================================
class Feeder : public threaded {
private:
  /// the task to perform
  const function<void()> task_;

public:
  Feeder(function<void()> task) : task_(std::move(task)) {};
protected:
  void worktask() const {task_();};
};

class Consumer : public threaded {
private:
  /// the task to perform
  function<void()> task_;

public:
  explicit Consumer(function<void()> task) : task_(std::move(task)) {};
protected:
  void worktask() const {(task_)();};
};



int main(int argc, char **argv) {
  auto machine = construct_algo(2., 0.5);

  LOG_INFO("Generate blibs");

  Generator gen(20., 10, 500);

  struct feeder_task {
    const Generator* gen_ptr_;
    TDBScan_AsyncMachine<ScalarBlib>* const machine_ptr_;
    feeder_task(const Generator* gen_ptr, TDBScan_AsyncMachine<ScalarBlib>* machine_ptr)
    : gen_ptr_(gen_ptr), machine_ptr_(machine_ptr) {};
    void operator()() const {machine_ptr_->FeedBlib(gen_ptr_->next());}
  } _feeder_task(&gen, &machine);


  // Feeder feeder(_feeder_task);
  //
  // auto consumer_task = [&machine]() {
  //   machine.ObtainCluster();
  // };
  //
  // Consumer consumerA(consumer_task), consumerB(consumer_task);

  //starting all parts of the machine
  machine.start();
  // consumerA.start();
  // consumerB.start();
  // feeder.start();

  std::this_thread::sleep_for(1s);
  machine.stop();

  //shutdown of all parts
}
