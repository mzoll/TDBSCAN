//
// Created by marcel on 16.04.21.
//

#include "external/common_clib/Semaphore.h"

#include <mutex>

using namespace std;

namespace common_clib::threadsafe {

Semaphore::Semaphore() noexcept
    : counter(0L),
      external_interrupt_set_counter_(0L),
      internal_interrupt_set_(false){};

Semaphore::~Semaphore() noexcept {
  std::lock_guard lock(mutex);
  cond.notify_all();
};

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


int Semaphore::consume_one() {
  if (is_blocked())
    throw interrupt_exception();
  std::unique_lock lock(mutex);
  lock.lock();
  if (counter.load() == 0) {
    cond.wait(lock, [&] {
      if (is_blocked())
        throw interrupt_exception();
      return counter > 0;
    });
  }
  --counter;
  return 1;
}


int Semaphore::consume_all() {
  if (is_blocked())
    throw interrupt_exception();
  std::lock_guard lock(mutex);
  auto rslt = counter.load();
  counter = 0;
  return rslt;
}


void Semaphore::set_interrupt() noexcept {
  std::lock_guard lock(mutex);
  internal_interrupt_set_.store(true);
  cond.notify_all();
};

void Semaphore::reset_interrupt() noexcept {
  std::lock_guard lock(mutex);
  internal_interrupt_set_.store(false);
};

void Semaphore::set_one_external_interrupt() noexcept {
  std::lock_guard lock(mutex);
  ++external_interrupt_set_counter_;
  cond.notify_all();
};

void Semaphore::reset_one_external_interrupt() noexcept {
  std::lock_guard lock(mutex);
  --external_interrupt_set_counter_;
};

}  // namespace common_clib::threadsafe
