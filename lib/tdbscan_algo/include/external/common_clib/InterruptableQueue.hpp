//
// Created by marcel on 16.04.21.
//
// inspiration from:
// https://stackoverflow.com/questions/27341029/what-is-the-best-way-to-wait-on-multiple-condition-variables-in-c11

/*
 * Defines an Interruptable Queue that can be used as an out-of-the-box solution
 * for access within a asynchronous context.
 * Guards and rails are put in place internally, so no need to manage critical sections of acces from the outside.
 */


#ifndef COMMONCLIB_INTERRUPTABLEQUEUE_H
#define COMMONCLIB_INTERRUPTABLEQUEUE_H

#include <functional>
#include <mutex>
#include <queue>

#include "external/common_clib/Semaphore.h"
#include "external/common_clib/trivial_logging.h"

namespace common_clib::threading {

/**
 * Asynchronous Queue class that works with a semaphore in order to allow an externally triggered interrupt-signal
 * to block. Lazy waits with an surfacing exception if either inlet or Outlet is closed
 * @tparam Tvalue the value-type held in the queue
 *
 * @example
 *   auto q = InterruptableQueue<int>()
 *
 *   std::thread pushing_thread([&q](){ std::thread::sleep_for(std::chrono::seconds(1)); q.push(42);})
 *   std::thread interrupting_thread([&q](){ std::thread::sleep_for(std::chrono::seconds(2)); q.block();})
 *
 *   q.push(10); // instantaneously put element '10' on the queue
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
  mutable common_clib::threadsafe::Semaphore semaphore_;  // <- coordinates 'post_one/consume_one' access by representing the queue's size

  /// represents the inlet and the outlet of the Queue
  class IOlet : public threadsafe::Interruptable {
    ///state if the iolet is closed or open for push or retrieval
    bool closed_{false};
  public:
    /// an Error that can be thrown when the Interrupt has been triggered
    struct IOletClosed : std::exception {
    private:
      const std::string message_;
    public:
      // Constructor accepting const char*
      IOletClosed(const char* msg) :
      message_(msg) {};

      [[nodiscard]] char const* what() const noexcept override {
        return message_.c_str();
      }
      ~IOletClosed() override = default;
    };

  public:
    void close(const std::function<void ()>* cb = nullptr) {
      closed_=true;
      set_internal_block();
      if (cb)
        (*cb)();
    };
    void open(const std::function<void ()>* cb=nullptr) {
      closed_=false;
      reset_internal_block();
      if (cb)
        (*cb)();
    };
    bool is_closed() const noexcept {return closed_;};
  } inlet_, outlet_;
public:
  using IOletClosed = IOlet::IOletClosed;

public:
  /// constructor
  InterruptableQueue() noexcept;
  /**
   * put a new entry on the queue; use move semantics
   * @note asynchronous call, not blocking
   * @exception IOletClosed thrown when the inlet has been closed
   */
  void push(value_type&& element);
  /**
   * put a new entry on the queue
   * @note asynchronous call, not blocking
   * @exception IOletClosed thrown when the inlet has been closed
   */
  void push(const value_type& element);
  /**
   * get a the oldest (first) entry from the queue and pop it, wait if not available intermediately
   * @note asynchronous call, blocking until an element is available
   * @exception IOletClosed thrown when the outlet has been closed
   * @throws common_clib::threadsafe::Interruptable::interrupt_set thrown when in wakeup fkt evaluates true
   * @return the oldest element of the queue
   * @example 
   *    bool end_execution{false};
   *    InterruptableQueue q;
   *    std::thread t([&]() {
   *      while (true) {
   *        try {
   *          q.pop([&](){return end_execution})
   *        } catch (interrupt_set) {
   *          //resurfacing from waiting 
   *          break;
   *        }   
   *      }
   *      // reaching here, means thread ends execution;
   *    });
   *    
   *    // t is thread; is blocked by waiting for consumption
   *    std::this_thread::sleep_for(std::chrono::seconds(1));
   *    end_execution = true; // will cause on next evaluation of _predicate_ in waiting to throw the `interrupt_set` 
   *    q.interrupt(); // triggers the wake-up from waiting and evaluation of _predicate_
   *    // sub-call `pop()` is interrupted and resurfaces to main thread control-loop, which is ending thread
   *    t.join();
   */
  [[nodiscard]] value_type pop(const std::function<bool()>& fkt = [](){return false;});
  /**
   * get all the entries that are currently on the queue
   * @exception IOletClosed thrown when the outlet has been closed
   * @return a vector of all remaining elements on the queue
   */
  [[nodiscard]] std::vector<value_type> exhaust();
  /// how many elements are on the ? Note that the state is highly volatile
  [[nodiscard]] size_t size() const noexcept;
  /// is the queue empty? Note that the state is highly volatile
  [[nodiscard]] bool empty() const noexcept;

  /// block the queue, all waiting calls will surface by throwing
  void close() noexcept;
  /// block the queue, all waiting calls will surface by throwing
  void close_inlet() noexcept;
  /// block the queue, all waiting calls will surface by throwing
  void close_outlet() noexcept;

  /** probe if queue has closed inlet or outlet
   * @return true if closed
   */
  bool is_closed() const noexcept;
  /// probe if queue's inlet is closed
  bool inlet_closed() const noexcept;
  /// probe if queue's outlet is closed
  bool outlet_closed() const noexcept;

  /// reopen any closed inlet or outlet
  void reopen() noexcept;
  /// reopen a closed inlet
  void reopen_inlet() noexcept;
  /// reopen a closed outlet
  void reopen_outlet() noexcept;

  /// purge all elements from the queue; this is done by blocking and unblocking the queue for purge
  void purge();

  /// get the block as an reference [that can be copied and proliferated]
  threadsafe::Interrupt<InterruptableQueue<Tvalue>>
  get_inlet_interrupt() const;
  /// get the block as an reference [that can be copied and proliferated]
  threadsafe::Interrupt<InterruptableQueue<Tvalue>>
  get_outlet_interrupt() const;

  /// trigger the interrupt on inlet and outlet
  void interrupt() const noexcept;

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
  if (inlet_.is_closed())
    throw IOletClosed("Inlet closed");
  std::lock_guard lock(push_mutex_);
  queue_.push(std::move(element));
  semaphore_.post_one();
}

template <class Tvalue>
void InterruptableQueue<Tvalue>::push(const value_type& element) {
  if (inlet_.is_closed())
    throw IOletClosed("Inlet closed");
  std::lock_guard lock(push_mutex_);
  queue_.push(element);
  semaphore_.post_one();
}

template <class Tvalue>
InterruptableQueue<Tvalue>::value_type
InterruptableQueue<Tvalue>::pop(const std::function<bool()>& fkt) {
  if (outlet_.is_closed())
    throw IOletClosed("Outlet is closed");
  try {
    semaphore_.consume_one(fkt);  // this will block until there is something to consume
  } catch (threadsafe::Interruptable::interrupt_set& e) {
    throw;
  } catch (threadsafe::Interruptable::block_set& e) {
    throw IOletClosed("Outlet has closed");
  }
  std::lock_guard pop_front_lock(pop_front_mutex_);
  auto element(std::move(queue_.front()));
  queue_.pop();
  return element;
}


template <class Tvalue>
std::vector<typename InterruptableQueue<Tvalue>::value_type>
InterruptableQueue<Tvalue>::exhaust() {
  std::vector<value_type> rslt_vec;
  if (outlet_.is_closed())
    throw IOletClosed("Outlet closed");
  unsigned int n_many;
  try {
    n_many = semaphore_.consume_all(); // blocking
  } catch (threadsafe::Interruptable::block_set& e) {
    throw IOletClosed("Outlet just closed");
  }

  // NOTE this might not work as expected, because consume and pop must be within one critical section, otherwise
  // another thread might steal the elements from the front of the queue, as new elements are pushed at the end
  std::lock_guard pop_front_lock(pop_front_mutex_);
  rslt_vec.reserve(n_many);
  for (unsigned long int i = 0; i < n_many; i++) {
    rslt_vec.emplace_back(std::move(queue_.front()));
    queue_.pop();
  }
  return rslt_vec;
}

template <class Tvalue>
void InterruptableQueue<Tvalue>::purge() {
  close();
  const auto n_many = semaphore_.consume_all();
  while (!queue_.empty()) {
    queue_.pop();
  }
  reopen();
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
void InterruptableQueue<Tvalue>::close() noexcept {
  // triggering the interrupt wakes up all waiting threads with an exception
  close_inlet();
  close_outlet();
}

template <class Tvalue>
void InterruptableQueue<Tvalue>::close_inlet() noexcept {
  // triggering the interrupt wakes up all waiting threads with an exception
  std::lock_guard lock(mutex_);
  inlet_.close();
}

template <class Tvalue>
void InterruptableQueue<Tvalue>::close_outlet() noexcept {
  // triggering the interrupt wakes up all waiting threads with an exception
  std::lock_guard lock(mutex_);
  outlet_.close();
  semaphore_.set_block();  //TODO, do we need to absorb this in into outlet?
}


template <class Tvalue>
bool InterruptableQueue<Tvalue>::is_closed() const noexcept {
  return inlet_closed() || outlet_closed();
}

template <class Tvalue>
bool InterruptableQueue<Tvalue>::inlet_closed() const noexcept {
  return inlet_.is_closed();
}

template <class Tvalue>
bool InterruptableQueue<Tvalue>::outlet_closed() const noexcept {
  return outlet_.is_closed();
}

template <class Tvalue>
void InterruptableQueue<Tvalue>::reopen() noexcept {
  inlet_.open();
  outlet_.open();
}

template <class Tvalue>
void InterruptableQueue<Tvalue>::reopen_inlet() noexcept {
  std::lock_guard lock(mutex_);
  inlet_.open();
}

template <class Tvalue>
void InterruptableQueue<Tvalue>::reopen_outlet() noexcept {
  std::lock_guard lock(mutex_);
  outlet_.open();
  semaphore_.reset_block();
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

template <class Tvalue>
void InterruptableQueue<Tvalue>::interrupt() const noexcept {
  semaphore_.interrupt();
}

}  // namespace common_clib::threading

#endif  // COMMONCLIB_INTERRUPTABLEQUEUE_H
