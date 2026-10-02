//
// Created by mzoll on 01.04.10.
//

#ifndef COMMONCLIB__CORE__STOPWATCH_H
#define COMMONCLIB__CORE__STOPWATCH_H

#include <chrono>
#include <sstream>
#include <string>
#include <vector>
#include <numeric>

// based on post of Francis Cugler on stackoverflow:
// https://stackoverflow.com/questions/22387586/measuring-execution-time-of-a-function-in-c

namespace common_clib {

/**
 * A simple Stopwatch that measures execution times.
 *
 * Embed this class in your code and use it in execution-time critical contexts. The stopwatch can either be started and stopped manually,
 * by holding onto a reference to the object, or it can be run at instantiation time with automatic stop and report at the end of scope,
 * which makes it very flexible in its use.
 *
 * @tparam Resolution the resolution requested, must derive from `std::chrono`-types, default `std::chrono::milliseconds`
 */
template <class Resolution = std::chrono::milliseconds>
class Stopwatch {
  

 public:
  enum class policy {
    defer = 0,
    start = 1
  };

 public:
  using timeunit_t = int;  // the type of `operator-(time - time).count()`
  using Clock =
      std::conditional_t<std::chrono::high_resolution_clock::is_steady,
                         std::chrono::high_resolution_clock,
                         std::chrono::steady_clock>;

 private:
  Clock::time_point mStart_ = Clock::now();
  Clock::time_point mLap_ = Clock::now();

  timeunit_t elapsed_total_{0};
  timeunit_t elapsed_lap_{0};
  std::vector<timeunit_t> laps_;

  bool started_ = false;
  bool paused_ = false;
  bool stopped_ = false;

  std::string name_;

 public:  //helper
  ///helper function to format some output
  [[nodiscard]] std::string timeunitString() const;

 public:  // Ctor/dtor
  /**
   * Construct an anonymous stopwatch with a policy.
   * @param ib the policy, either 'start' (default) or 'deferred'. In case of start the clock, starts ticking immediately,
   * in case of 'deferred' the clock needs to be started manually by calling `start()`.
   */
  explicit Stopwatch(policy ib = policy::start);
  /**
   * Construct a named stopwatch with a policy.
   * @param name a name to the StopWatch; Will be used for reporting
   * @param ib the policy, either 'start' (default) or 'deferred'. In case of start the clock, starts ticking immediately,
   * in case of 'deferred' the clock needs to be started manually by calling `start()`.
   */
  explicit Stopwatch(std::string name, policy ib = policy::start);

 public:  //methods
  /** start the stopwatch.
   *
   * @throws runtime_error in case the stopwatch has already been started; use `restart()` for nothrow.
   */
  void start();
  
  /** start the stopwatch; allows restarting
   *
   * @throws runtime_error in case the stopwatch has already been started; use `restart()` for nothrow.
   */
  void pause();

  /** reset the stopwatch. must be called after `stop()` for (re)starting.
   *
   * * @note does not throw, even if the stopwatch has not been started before.
   */
  void reset();

  /** restart the stopwatch
   *
   * @note does not throw, even if the stopwatch has not been started before
   */
  void restart();

  ///obtain the lap-time, aka the time since the last time the lap funktion had been used or since start
  [[maybe_unused]] timeunit_t lap();
  ///obtain the time since start and stop the stopwatch and
  [[maybe_unused]] timeunit_t stop();

  ///obtain the elapsed time since start (without stopping the stopwatch)
  [[nodiscard]] timeunit_t time();

  /**
   * take the lap time and make a formatted printout of the elapsed time.
   *
   ** @throws runtime_error in case the stopwatch was not started before calling this function
   */
  [[nodiscard]] std::string lap_report_elapsed(const std::string& tp_message);
  /**
   * stop the clock and make a formatted printout of the elapsed time since start.
   *
   * @param tp_message a custom message attached to this print out
   * @throws runtime_error in case the stopwatch was not stopped before calling this function
   */
  [[nodiscard]] std::string stop_report_elapsed(const std::string& tp_message);

  /** make a formatted printout of the elapsed time; does NOT stop the clock.
   *
   * @param tp_message a custom message attached to this print out
   */
  [[nodiscard]] std::string report_elapsed(const std::string& tp_message);

  /** make a formatted printout of the elapsed time; does NOT stop the clock.
   *
   * @param tp_message a custom message attached to this print out
   */
  [[nodiscard]] std::string report_avglap(const std::string& tp_message);
};

}  // namespace common_clib

//=========================
//====== DEFINITIONS ======
//=========================

namespace common_clib {

// ExecutionTimer() = default;
template <class Resolution>
Stopwatch<Resolution>::Stopwatch(const policy ib)
    : Stopwatch("", ib) {};

template <class Resolution>
Stopwatch<Resolution>::Stopwatch(std::string name, const policy ib )
    : name_(name) {
  if (ib == policy::start) {
    mStart_ = Clock::now();
    mLap_ = mStart_;
    started_ = true;
  }
};

template <class Resolution>
void Stopwatch<Resolution>::reset() {
  mStart_ = Clock::now();
  mLap_ = mStart_;
  elapsed_total_ = 0;
  elapsed_lap_ = 0;
  laps_.clear();
  started_ = false;
}

template <class Resolution>
void Stopwatch<Resolution>::start() {
  if (started_)
    throw std::runtime_error("Stopwatch already started");

  mStart_ = Clock::now();
  mLap_ = mStart_;
  started_ = true;
}

template <class Resolution>
void Stopwatch<Resolution>::pause() {
  if (! started_)
    throw std::runtime_error("Cannot pause the Stopwatch before it has been started");
  const auto end_time = Clock::now();

  elapsed_total_ += std::chrono::duration_cast<Resolution>(end_time - mStart_).count();
  elapsed_lap_ += std::chrono::duration_cast<Resolution>(end_time - mStart_).count();

  paused_ = true;
}


template <class Resolution>
void Stopwatch<Resolution>::restart() {
  if (! stopped_)
    throw std::runtime_error("Cannot restart Stopwatch after having stopped; should have used `paused()` instead");

  mStart_ = Clock::now();
  mLap_ = mStart_;
  started_ = true;
  paused_ = false;
}

template <class Resolution>
Stopwatch<Resolution>::timeunit_t Stopwatch<Resolution>::lap() {
  if (! started_)
    throw std::runtime_error("Stopwatch was never started");
  if (stopped_)
    return elapsed_lap_;
  const auto laptstart_timestamp = mLap_;
  mLap_ = Clock::now();
  const auto this_lap = std::chrono::duration_cast<Resolution>(mLap_ - laptstart_timestamp).count() + elapsed_lap_;

  laps_.push_back(elapsed_lap_);
  elapsed_lap_ = 0;

  return this_lap;
}


template <class Resolution>
Stopwatch<Resolution>::timeunit_t Stopwatch<Resolution>::stop() {
  if (! started_)
    throw std::runtime_error("Cannot stop the Stopwatch before it has been started");
  if (stopped_ or paused_) {
    stopped_ = true;
    return elapsed_total_;
  }

  const auto end_time = Clock::now();
  stopped_ = true;

  elapsed_total_ += std::chrono::duration_cast<Resolution>(end_time - mStart_).count();
  elapsed_lap_ += std::chrono::duration_cast<Resolution>(end_time - mStart_).count();

  laps_.push_back(elapsed_lap_);

  return elapsed_total_;
}

template <class Resolution>
Stopwatch<Resolution>::timeunit_t Stopwatch<Resolution>::time() {
  if (! started_)
    throw std::runtime_error("Cannot obtain a time reading before the Stopwatch has been started");
  if (stopped_)
    return elapsed_total_;

  const auto current_time = Clock::now();
  return std::chrono::duration_cast<Resolution>(current_time - mStart_).count() + elapsed_total_;
}


template <class Resolution>
std::string Stopwatch<Resolution>::lap_report_elapsed(const std::string& tp_message) {
  const auto _lap = lap();
  std::ostringstream ss;
  ss << "Stopwatch " << name_ << ": " << tp_message << " :: " << _lap << " " << timeunitString();
  return ss.str();
}


template <class Resolution>
std::string Stopwatch<Resolution>::stop_report_elapsed(const std::string& tp_message) {
  const auto _stop = stop();
  std::ostringstream ss;
  ss << "Stopwatch " << name_ << ": " << tp_message << " :: " << _stop << " " << timeunitString();
  return ss.str();
}

template <class Resolution>
std::string Stopwatch<Resolution>::report_elapsed(const std::string& tp_message) {
  const auto _time =  time();
  std::ostringstream ss;
  ss << "Stopwatch " << name_ << ": " << tp_message << " :: " << _time << " " << timeunitString();
  return ss.str();
}

template <class Resolution>
std::string Stopwatch<Resolution>::report_avglap(const std::string& tp_message) {

  const auto _sum = std::accumulate(std::begin(laps_), std::end(laps_), 0.0);
  const double avglap = _sum / laps_.size();

  std::ostringstream ss;
  ss << "Stopwatch " << name_ << ": " << tp_message << " :: " << avglap << " " << timeunitString();
  return ss.str();
}


template <class Resolution>
std::string Stopwatch<Resolution>::timeunitString() const {
  if (std::is_same<Resolution, std::chrono::seconds>::value) {
    return "s";
  }
  if (std::is_same<Resolution, std::chrono::milliseconds>::value) {
    return "ms";
  }
  if (std::is_same<Resolution, std::chrono::microseconds>::value) {
    return "us";
  }
  if (std::is_same<Resolution, std::chrono::nanoseconds>::value) {
    return "ns";
  }

  return " (unit unknown)";
}

}  // namespace common_clib

#endif // COMMONCLIB__CORE__STOPWATCH_H
