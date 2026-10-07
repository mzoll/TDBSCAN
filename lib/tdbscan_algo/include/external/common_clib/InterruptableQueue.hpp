//
// Created by marcel on 16.04.21.
//
// inspiration from:
// https://stackoverflow.com/questions/27341029/what-is-the-best-way-to-wait-on-multiple-condition-variables-in-c11

/*
 * Defines an Interruptable Queue that can be used as an out of the box solution for access within a asynchronous context.
 * Guards and rails are pu in place internally, so no need to manage this outside
 *
 */


#ifndef COMMONCLIB_INTERRUPTABLEQUEUE_H
#define COMMONCLIB_INTERRUPTABLEQUEUE_H

#include <functional>
#include <mutex>
#include <queue>

#include "external/common_clib/Semaphore.h"

namespace common_clib::threading {

/// an Error to be thrown when an item should have been popped from the queue, but it was empty
struct QueueEmptyError : std::exception {
  [[nodiscard]] char const* what() const noexcept override {
    return "The Queue was queried for an item, when there wasn't one available";
  }
  ~QueueEmptyError() override = default;
};

/**
 * Asynchronous Queue class that works with a semaphore in order to allow an externally
 triggered interrupt-signal to block
 * lazy waits with an surfacing exception
 * @tparam Tvalue the value-type held in the queue
 *
 * @example
 *   auto q = InterruptableQueue<int>()
 *
 *   std::thread pushing_thread([&q](){ std::thread::sleep_for(std::chrono::seconds(1)); q.push(42);})
 *   std::thread interrupting_thread([&q](){ std::thread::sleep_for(std::chrono::seconds(2)); q.block();})
 *
 *   q.push(10); // instantaniously put element '10' on the queue
 *   auto element0 = q.pop(); // <- =10; [more or less instantaneous]
 *   auto element1 = q.pop_or_wait(); // <- =42; blocking until element becomes available aka. `pushing_thread` has pushed
 *   try
 *      auto element2 = q.pop_or_wait(); // blocking, as now element is available on queue; after 2 seconds: resurfaces with InterruptIsSetError, because interrupting_thread's action
 *   catch (InterruptIsSetError) {
 *     pass;
 *   }
 *   q.release_block()
 */
template <class Tvalue>
class InterruptableQueue {
 public:
  /// the type of the stored elements
  using value_type = Tvalue;

  /// an Error that can be thrown when the Interrupt has been triggered
  struct InterruptIsSet : std::exception {
    [[nodiscard]] char const* what() const noexcept override {
      return "The Queue was interrupted";
    }
    ~InterruptIsSet() override = default;
  };

 private:
  /// the most elemental queue-container
  std::queue<value_type> queue_;
  /// coordinates simultaneous access between threads, when in serial access is required
  std::mutex mutex_;

  ///needed to lock access on the queue at the from together
  std::mutex push_mutex_;
  ///needed to lock access on the queue at the from together
  std::mutex pop_front_mutex_;

  /// coordinates access external interrupts
  std::mutex io_interrupt_mutex_;

  /// an atomic counter, that more or less represents the size of the queue and
  /// lets poping-threads wait
  mutable common_clib::threadsafe::Semaphore semaphore_;  // <- coordinates access by representing the queue's size

  class _IOlet : public threadsafe::Interruptable {
    bool closed_{false};
  public:
    void block(const std::function<void ()>* cb = nullptr) {
      closed_=true;
      set_internal_interrupt();
      if (cb)
        (*cb)();
    };
    void release(const std::function<void ()>* cb=nullptr) {
      closed_=false;
      reset_internal_interrupt();
      if (cb)
        (*cb)();
    };
    bool is_blocked() const noexcept {return closed_;};
  } inlet_, outlet_;


 public:
  /// constructor
  InterruptableQueue() noexcept;
  /// put a new entry on the queue; use move semantics
  void push(value_type&& element);
  /// put a new entry on the queue
  void push(const value_type& element);

  /**
   * get a the oldest (first) entry from the queue and pop it, wait if not available intermediately
   * @return the oldest element of the queue
   */
  [[nodiscard]] value_type pop();
  /// get all the entries that are currently on the queue
  [[nodiscard]] std::vector<value_type> exhaust();
  /// how many elements are on the ? Note that the state is highly volatile
  [[nodiscard]] size_t size() const noexcept;
  /// is the queue empty? Note that the state is highly volatile
  [[nodiscard]] bool empty() const noexcept;

  /// block the queue, all waiting calls will surface by throwing
  void block() noexcept;
  /// block the queue, all waiting calls will surface by throwing
  void block_inlet() noexcept;
  /// block the queue, all waiting calls will surface by throwing
  void block_outlet() noexcept;

  /// probe if queue is blocked
  bool is_blocked() const noexcept;
  /// probe if queue is blocked
  bool inlet_blocked() const noexcept;
  /// probe if queue is blocked
  bool outlet_blocked() const noexcept;

  /// reset the block
  void release_block() noexcept;
/// reset the block
  void release_inlet_block() noexcept;
/// reset the block
  void release_outlet_block() noexcept;

  /// purge all elements from the queue; this is done by blocking and unblocking the queue for purge
  void purge();

  /// get the block as an reference [that can be copied and proliferated]
  threadsafe::Interrupt<InterruptableQueue<Tvalue>>
  get_inlet_interrupt() const;
  /// get the block as an reference [that can be copied and proliferated]
  threadsafe::Interrupt<InterruptableQueue<Tvalue>>
  get_outlet_interrupt() const;

  friend common_clib::threadsafe::Interrupt<InterruptableQueue>;
};

}  // namespace common_clib::threading

// ============================
// =====  IMPLEMENTATIONS =====
// ============================

namespace common_clib::threading {

template <class Tvalue>
InterruptableQueue<Tvalue>::InterruptableQueue() noexcept : queue_() {}

template <class Tvalue>
void InterruptableQueue<Tvalue>::push(value_type&& element) {
  if (inlet_.is_blocked())
    throw InterruptIsSet();
  std::lock_guard lock(push_mutex_);
  queue_.push(std::move(element));
  semaphore_.post_one();
}

template <class Tvalue>
void InterruptableQueue<Tvalue>::push(const value_type& element) {
  if (inlet_.is_blocked())
    throw InterruptIsSet();
  std::lock_guard lock(push_mutex_);
  queue_.push(element);
  semaphore_.post_one();
}

template <class Tvalue>
InterruptableQueue<Tvalue>::value_type
InterruptableQueue<Tvalue>::pop() {
  if (outlet_.is_blocked())
    throw InterruptIsSet();
  try {
    semaphore_.consume_one();  // this will block if there is currently nothing to consume

  } catch (threadsafe::interrupt_exception& e) {
    throw InterruptIsSet();
  } catch (...) {
    throw InterruptIsSet();
  }

  std::lock_guard pop_front_lock(pop_front_mutex_);
  auto element(std::move(queue_.front()));
  queue_.pop();
  return element;
}


template <class Tvalue>
std::vector<typename InterruptableQueue<Tvalue>::value_type>
InterruptableQueue<Tvalue>::exhaust() {
  if (outlet_.is_blocked())
    throw InterruptIsSet();
  unsigned int n_many;
  try {
    n_many = semaphore_.consume_all(); // blocking
  } catch (threadsafe::interrupt_exception& e) {
    throw InterruptIsSet();
  }

  std::vector<value_type> rslt_vec;
  rslt_vec.reserve(n_many);
  std::lock_guard pop_front_lock(pop_front_mutex_);
  for (unsigned long int i = 0; i < n_many; i++) {
    rslt_vec.emplace_back(std::move(queue_.front()));
    queue_.pop();
  }
  return rslt_vec;
}

template <class Tvalue>
void InterruptableQueue<Tvalue>::purge() {
  block();
  const auto n_many = semaphore_.consume_all();
  while (!queue_.empty()) {
    queue_.pop();
  }
  release_block();
}


template <class Tvalue>
size_t InterruptableQueue<Tvalue>::size() const noexcept {
  return semaphore_.load();
}

template <class Tvalue>
bool InterruptableQueue<Tvalue>::empty() const noexcept {
  return semaphore_.load() == 0;
}

template <class Tvalue>
void InterruptableQueue<Tvalue>::block() noexcept {
  // triggering the interrupt wakes up all waiting threads with an exception
  block_inlet();
  block_outlet();
}

template <class Tvalue>
void InterruptableQueue<Tvalue>::block_inlet() noexcept {
  // triggering the interrupt wakes up all waiting threads with an exception
  std::lock_guard lock(mutex_);
  inlet_.block();
}

template <class Tvalue>
void InterruptableQueue<Tvalue>::block_outlet() noexcept {
  // triggering the interrupt wakes up all waiting threads with an exception
  std::lock_guard lock(mutex_);
  outlet_.block();
  semaphore_.set_interrupt();  //TODO, do we need to absorb this in into outlet?
}


template <class Tvalue>
bool InterruptableQueue<Tvalue>::is_blocked() const noexcept {
  return inlet_blocked() || outlet_blocked();
}

template <class Tvalue>
bool InterruptableQueue<Tvalue>::inlet_blocked() const noexcept {
  return inlet_.is_blocked();
}

template <class Tvalue>
bool InterruptableQueue<Tvalue>::outlet_blocked() const noexcept {
  return outlet_.is_blocked();
}

template <class Tvalue>
void InterruptableQueue<Tvalue>::release_block() noexcept {
  inlet_.release();
  outlet_.release();
}

template <class Tvalue>
void InterruptableQueue<Tvalue>::release_inlet_block() noexcept {
  std::lock_guard lock(mutex_);
  inlet_.release();
}

template <class Tvalue>
void InterruptableQueue<Tvalue>::release_outlet_block() noexcept {
  std::lock_guard lock(mutex_);
  outlet_.release();
  semaphore_.reset_interrupt();
}

/// get the block as an reference [that can be copied and proliferated]
template <class Tvalue>
common_clib::threadsafe::Interrupt<InterruptableQueue<Tvalue>>
InterruptableQueue<Tvalue>::get_inlet_interrupt() const {
  return inlet_.get_interrupt();
}

/// get the block as an reference [that can be copied and proliferated]
template <class Tvalue>
common_clib::threadsafe::Interrupt<InterruptableQueue<Tvalue>>
InterruptableQueue<Tvalue>::get_outlet_interrupt() const {
  return outlet_.get_interrupt();
}



}  // namespace common_clib::threading

#endif  // COMMONCLIB_INTERRUPTABLEQUEUE_H
