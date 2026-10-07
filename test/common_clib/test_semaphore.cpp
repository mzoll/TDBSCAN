//
// Created by marcel on 23.02.21.
//


#include <gtest/gtest.h>

#include "external/common_clib/Semaphore.h"

#include <chrono>
#include <exception>
#include <future>
#include <thread>


using namespace std;
using namespace common_clib;
using namespace common_clib::threadsafe;

TEST(Semaphore, serial_access) {
  Semaphore s;

  EXPECT_FALSE(s.is_blocked());

  EXPECT_EQ(s.load(), 0);

  EXPECT_NO_THROW(s.post_one());
  EXPECT_EQ(s.load(), 1);

  // set the interrupt
  EXPECT_NO_THROW(s.set_interrupt());
  EXPECT_TRUE(s.is_blocked());
  EXPECT_NO_THROW(s.post_one());
  EXPECT_EQ(s.load(), 2);

  // reset the interrupt
  EXPECT_NO_THROW(s.reset_interrupt());
  EXPECT_FALSE(s.is_blocked());
  EXPECT_NO_THROW(s.post_one());
  EXPECT_EQ(s.load(), 3);

  // consume from here
  EXPECT_NO_THROW(s.consume_one());
  EXPECT_EQ(s.load(), 2);

  // set the interrupt
  EXPECT_NO_THROW(s.set_interrupt());
  EXPECT_TRUE(s.is_blocked());
  EXPECT_ANY_THROW(s.consume_one());
  EXPECT_EQ(s.load(), 2);

  // reset the interrupt
  EXPECT_NO_THROW(s.reset_interrupt());
  EXPECT_FALSE(s.is_blocked());
  EXPECT_NO_THROW(s.consume_one());
  EXPECT_EQ(s.load(), 1);
}

TEST(Semaphore, parallel_access) {
  Semaphore s;

  std::promise<bool> p0;

  const auto task0 = [&s, &p0]() {
    try {
      s.consume_one();
    } catch (const interrupt_exception& e) {
      p0.set_exception(std::current_exception());
      return;
    }
    p0.set_value_at_thread_exit(true);
  };

  std::thread t0(task0);
  std::this_thread::sleep_for(std::chrono::milliseconds(500));
  auto f0 = p0.get_future();

  f0.wait_for(std::chrono::milliseconds(100));

  s.post_one();
  t0.join();
  EXPECT_TRUE(f0.get() == true);
  EXPECT_TRUE(s.load() == 0);

  // now lets, see if we can set an interrupt on the waiting
  std::promise<bool> p1;
  const auto task1 = [&s, &p1]() {
    try {
      s.consume_one();
    } catch (const interrupt_exception& e) {
      p1.set_exception(std::current_exception());
      return;
    }
    p1.set_value_at_thread_exit(true);
  };

  std::thread t1(task1);
  std::this_thread::sleep_for(std::chrono::milliseconds(500));
  auto f1 = p1.get_future();

  s.set_interrupt();
  s.post_one();
  t1.join();
  EXPECT_ANY_THROW(f1.get());
  EXPECT_TRUE(s.load() == 1);
}

TEST(Semaphore, interrupts) {
  Semaphore s;

  // interrupts that might be distributed among esternal (otherwise inaccesible)
  //  tasks
  auto i0 = s.get_interrupt();
  auto i1 = s.get_interrupt();

  EXPECT_FALSE(s.is_blocked());

  // set the trigger on i0
  EXPECT_FALSE(i0.is_triggered());
  i0.trigger();
  EXPECT_TRUE(i0.is_triggered());
  EXPECT_TRUE(s.is_blocked());
  i0.trigger();  // triggering an already triggered Interrupt is inert

  // reset the i0 trigger
  i0.reset();
  EXPECT_FALSE(i0.is_triggered());
  EXPECT_FALSE(s.is_blocked());

  i0.reset();  // resetting an already reset Interrupt is inert

  // trigger interrupt from multiple sites
  EXPECT_FALSE(s.is_blocked());
  i0.trigger();
  i1.trigger();
  EXPECT_TRUE(s.is_blocked());
  i0.reset();
  EXPECT_TRUE(s.is_blocked());
  i1.reset();
  EXPECT_FALSE(s.is_blocked());
}