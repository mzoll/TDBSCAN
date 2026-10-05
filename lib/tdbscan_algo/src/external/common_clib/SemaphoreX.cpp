//
// Created by marcel on 16.04.21.
//

#include "external/common_clib/SemaphoreX.h"

#include <mutex>

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


int Semaphore::consume_one() {
  if (is_blocked())
    throw interrupt_exception();
  std::unique_lock lock(mutex);

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

}  // namespace common_clib::threadsafe
