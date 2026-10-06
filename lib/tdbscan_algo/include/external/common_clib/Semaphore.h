//
// Created by marcel on 16.04.21.
//
// A class that allows to have multiple consumers wait,
// until a notification signal is send; very similar to the mechanics of
// conditional variable. However this construct allows the obtaining of a
// interrupt signal, which when _trigggered_ will resurface on the wait with an
// _interrupt_exception_.

#ifndef COMMONCLIB__THREADSAFE__SEMAPHORE_H
#define COMMONCLIB__THREADSAFE__SEMAPHORE_H

#include <atomic>
#include <condition_variable>
#include <mutex>

#include "interrupt.h"

namespace common_clib::threadsafe {

/**
 * A class that allows to have multiple consumers waiting,
 * until a notification signal is send; very similar to the mechanics of
 * conditional variable. However this construct allows the obtaining of a block
 * signal, which when _triggered_ will resurface on the wait with an
 * _interrupt_exception_.
 *
 * @example
 *    //posting and consuming
 *    using namespace std::chrono_literals;
 *
 *    Semaphore s;
 *
 *    std::thread consumer( [&s],{
 *      unsigned int counter = 0;
*        while (true) {
 *        try
   *        s.consume_one();
   *        ++counter;
   *      } catch (interrupt_exception) {
   *        std::cout << "consumption interrupted after " << counter << std::endl;
 *      }
 *    });
 *
 *    std::thread t0( [&s],{
 *      std::this_thread::sleep_for(1s);
 *      s.post_one();
 *    });
 *
*     std::thread t0( [&s],{
 *      std::this_thread::sleep_for(1s);
 *      s.post_many(2);
 *    });
 *
 *    std::this_thread::sleep_for(2s);
 *    s.set_interrupt(); // <<< consumption interrupted after 3
 *
 *    t0.join();
 *    consumer.join();
 *
 * @example
 *    // interrupt dangling thread
*     Semaphore s;
 *
 *    std::thread consumer( [&s]() {
*      unsigned int counter = 0;
*        while (true) {
 *        try
   *        s.consume_one();
   *        ++counter;
   *      } catch (interrupt_exception) {
   *        std::cout << "consumption interrupted after " << counter << std::endl;
 *      }
 *    });
 *
 *    auto i = s.get_interrupt();
 *
 *    // pass on the interrupt handle rather than to pass on the semaphore itself
 *    std::thread breaker([&i]() {
 *      std::this_thread::sleep_for(1s);
 *      i.trigger();
 *    });
 *
 *    consumer.join();
 *    breaker.join();
 */
class Semaphore : protected Interruptable {
 private:
  mutable std::mutex mutex;
  // the central counting unit
  std::atomic<long> counter;

 public:  // ctor
  Semaphore() noexcept;

 public:
  /**
   * post *one* element, allowing the consumption of *One* element; this is semi lock-free
   * [raise the counter by one and notify one waiting party (release the block on wait())]
   * @throws interrupt_exception if the semaphore has been blocked
   */
  void post_one();

  /**
   * post that *many* elements at once;  allowing the consumption of that *many* elements; this is semi lock-free
   * @param n_many post that many
   * @throws interrupt_exception if the semaphore has been blocked
   */
  void post_many(int n_many);

  /** consume *One* element if immediately available, or wait [block] until it becomes available.
   * @throws interrupt_exception if the semaphore has been blocked
   */
  [[maybe_unused]] int consume_one();

  /** consume *All* elements possible
   * @throws interrupt_exception if the semaphore has been blocked
   * @return the number of consumed elements
   */
  [[maybe_unused]] int consume_all();

  // convenience access
  /// get the number of consumable elements [volatile]
  unsigned long load() const { return counter.load(); }

  /// get the current status
  bool is_blocked() const {
    return is_interrupted();
  }

  /// toggle to blocked state, raising an exception on all waiting parties; no
  /// further consume operations are allowed until the reset_interrupt is
  /// called
  void set_interrupt() noexcept;

  /// toggle from blocked state to free running,
  void reset_interrupt() noexcept;

  // friend
  friend Interrupt<Semaphore>;
};

}  // namespace common_clib::threadsafe

#endif  // COMMONCLIB__THREADSAFE__SEMAPHORE_H
