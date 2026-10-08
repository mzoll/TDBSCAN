//
// Created by netsu on 08/10/2026.
//

#ifndef TDBSCAN__EXASYNC__AGENTS_H
#define TDBSCAN__EXASYNC__AGENTS_H

#include <functional>
#include <thread>
#include <utility>

namespace exasync {

using namespace std;
using namespace exasync;
using namespace tdbscan;


/**
 * A prototype for threaded things
 */
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

} // namespace exasync

#endif //TDBSCAN__EXASYNC__AGENTS_H
