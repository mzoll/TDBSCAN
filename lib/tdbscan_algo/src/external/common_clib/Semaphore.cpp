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

int Semaphore::consume_one(const function<bool()> *fkt) {
  if (is_blocked())
    throw block_set("consumption is interrupted");
  std::unique_lock lock(mutex);
  if (counter.load() == 0) {
    cond.wait(lock, [&] {
      if (fkt && (*fkt)())
        throw interrupt_set("interrupt was set");
      if (is_blocked())
        throw block_set("consumption is interrupted");
      return counter > 0;
    });
  }
  --counter;
  return 1;
}


int Semaphore::consume_all() {
  if (is_blocked())
    throw block_set("consumption is interrupted");
  std::lock_guard lock(mutex); //this is blocking until lock can be acquired
  if (is_blocked())
    throw block_set("consumption is interrupted");
  auto rslt = counter.load();
  counter = 0;
  return rslt;
}

unsigned long Semaphore::load() const noexcept {
  return counter.load();
}

void Semaphore::set_block() noexcept {
  set_internal_block();
};

void Semaphore::reset_block() noexcept {
  reset_internal_block();
};

void Semaphore::interrupt() noexcept {
  Interruptable::interrupt();
};



}  // namespace common_clib::threadsafe
