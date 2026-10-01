//
// Created by mzoll on 01.04.10.
//

#ifndef COMMONCLIB__CORE__STOPWATCH_H
#define COMMONCLIB__CORE__STOPWATCH_H

#include <chrono>
#include <sstream>
#include <string>

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
  using Clock =
      std::conditional_t<std::chrono::high_resolution_clock::is_steady,
                         std::chrono::high_resolution_clock,
                         std::chrono::steady_clock>;

 private:
  Clock::time_point mStart_ = Clock::now();
  Clock::time_point mLap_ = Clock::now();
  std::string name_;
  bool started_ = false;
 private:  //helper
  ///helper function to format some output
  [[nodiscard]] std::string PrintTimeunit() const;

 public:  //Ctor

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
  /// start the stopwatch if it has not already been started
  void start();
  /// restart the stopwatch
  void restart();
  ///obtain the lap-time, aka the time since the last time the lap funktion had been used or since start
  int lap();
  ///make a formatted printout of the elapsed time for this lap
  std::string lap_report(const std::string& tp_message);
  ///obtain the time since start and stop the stopwatch and
  int stop();
  ///make a formatted printout of the elapsed time since start
  std::string stop_report(const std::string& tp_message);
  ///obtain the elapsed time since start (without stopping the stopwatch)
  int time();
  ///make a formatted printout of the elapsed time
  std::string time_report(const std::string& tp_message);
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
void Stopwatch<Resolution>::start() {
  if (started_)
    throw std::runtime_error("Stopwatch already started");

  mStart_ = Clock::now();
  mLap_ = mLap_;
  started_ = true;
}

template <class Resolution>
void Stopwatch<Resolution>::restart() {
  mStart_ = Clock::now();
  mLap_ = mStart_;
  started_ = true;
}

template <class Resolution>
int Stopwatch<Resolution>::lap() {
  if (! started_)
    throw std::runtime_error("Cannot stop the Stopwatch before it has been started");
  const auto end_time = Clock::now();
  const auto start_time = mLap_;
  mLap_ = Clock::now();
  return std::chrono::duration_cast<Resolution>(end_time - start_time).count();

}

template <class Resolution>
std::string Stopwatch<Resolution>::lap_report(const std::string& tp_message) {
  if (! started_)
    throw std::runtime_error("Cannot stop the Stopwatch before it has been started");
  const auto end = lap();
  std::ostringstream ss;
  ss << "Stopwatch " << name_ << ": " << tp_message << " :: "
     << end << " " << PrintTimeunit();
  return ss.str();
}

template <class Resolution>
int Stopwatch<Resolution>::stop() {
  if (! started_)
    throw std::runtime_error("Cannot stop the Stopwatch before it has been started");
  const auto end_time = Clock::now();
  started_ = false;
  return std::chrono::duration_cast<Resolution>(end_time - mStart_).count();
}

template <class Resolution>
std::string Stopwatch<Resolution>::stop_report(const std::string& tp_message) {
  if (! started_)
    throw std::runtime_error("Cannot stop the Stopwatch before it has been started");

  const auto end = stop();
  std::ostringstream ss;
  ss << "Stopwatch " << name_ << ": " << tp_message << " :: "
     << end << " " << PrintTimeunit();
  return ss.str();
}

template <class Resolution>
int Stopwatch<Resolution>::time() {
  if (! started_)
    throw std::runtime_error("Cannot obtain a time reading before the Stopwatch has been started");
  const auto end_time = Clock::now();
  return std::chrono::duration_cast<Resolution>(end_time - mStart_).count();
}

template <class Resolution>
std::string Stopwatch<Resolution>::time_report(const std::string& tp_message) {
  if (! started_)
    throw std::runtime_error("Cannot obtain a time reading before the Stopwatch has been started");

  const auto end = time();
  std::ostringstream ss;
  ss << "Stopwatch " << name_ << ": " << tp_message << " :: "
     << end << " " << PrintTimeunit();
  return ss.str();
}


template <class Resolution>
std::string Stopwatch<Resolution>::PrintTimeunit() const {
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
