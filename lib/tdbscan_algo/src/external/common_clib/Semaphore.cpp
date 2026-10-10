//
// Created by marcel on 16.04.21.
//

#include "external/common_clib/Semaphore.h"

#include <mutex>
#include <functional>

#include "external/common_clib/trivial_logging.h"

using namespace std;

namespace common_clib::threadsafe {

Semaphore::Semaphore() noexcept
    : counter(0L) {};

void Semaphore::post_one() {
  std::lock_guard lock(mutex);
  ++counter;
  cond.notify_one();  // never throws
}

void Semaphore::post_many(const int n_many) {
  std::lock_guard lock(mutex);
  for (int n=0; n<n_many; n++) {
    ++counter;
    cond.notify_one();
  }
}

int Semaphore::consume_one(const function<void()> *fkt) {
  if (is_blocked())
    throw interrupt_exception("consumption is interrupted");
  std::unique_lock lock(mutex);
  if (counter.load() == 0) {
    LOG_FATAL("WAS here");
    cond.wait(lock, [&] {
      if (is_blocked())
        throw interrupt_exception("consumption is interrupted");
      if (fkt)
        (*fkt)();
      return counter > 0;
    });
  }
  --counter;
  return 1;
}


int Semaphore::consume_all() {
  if (is_blocked())
    throw interrupt_exception("consumption is interrupted");
  std::lock_guard lock(mutex);
  auto rslt = counter.load();
  counter = 0;
  return rslt;
}

unsigned long Semaphore::load() const noexcept {
  return counter.load();
}

void Semaphore::set_interrupt() noexcept {
  set_internal_block();
};

void Semaphore::reset_interrupt() noexcept {
  reset_internal_block();
};



}  // namespace common_clib::threadsafe
