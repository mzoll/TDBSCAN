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
  EXPECT_NO_THROW(s.set_block());
  EXPECT_TRUE(s.is_blocked());
  EXPECT_NO_THROW(s.post_one());
  EXPECT_EQ(s.load(), 2);

  // reset the interrupt
  EXPECT_NO_THROW(s.reset_block());
  EXPECT_FALSE(s.is_blocked());
  EXPECT_NO_THROW(s.post_one());
  EXPECT_EQ(s.load(), 3);

  // consume from here
  EXPECT_NO_THROW(s.consume_one({}));
  EXPECT_EQ(s.load(), 2);

  // set the interrupt
  EXPECT_NO_THROW(s.set_block());
  EXPECT_TRUE(s.is_blocked());
  EXPECT_ANY_THROW(s.consume_one({}));
  EXPECT_EQ(s.load(), 2);

  // reset the interrupt
  EXPECT_NO_THROW(s.reset_block());
  EXPECT_FALSE(s.is_blocked());
  EXPECT_NO_THROW(s.consume_one({}));
  EXPECT_EQ(s.load(), 1);
}

TEST(Semaphore, async_access) {
  //TODO find out if promise can be recycled

  Semaphore s;
  {
    std::promise<bool> p0;

    const auto task0 = [&s, &p0]() {
      try {
        s.consume_one({});  // blocks if there is nothing to consume
      } catch (const Interruptable::interrupt_set& e) {
        p0.set_exception(std::current_exception());
        return;
      } catch (const Interruptable::block_set& e) {
        p0.set_exception(std::current_exception());
        return;
      }
      p0.set_value_at_thread_exit(true);
    };

    std::thread t0(task0);  //this is executes but, as there is nothing to consume, it will hang
    std::this_thread::sleep_for(std::chrono::milliseconds(50));
    auto f0 = p0.get_future();

    f0.wait_for(std::chrono::milliseconds(10));  // waits on the arrival of the result, but it wont come, unless...

    s.post_one(); // something is actually pushed, waking up t0, and consumption continues
    t0.join();
    EXPECT_EQ(f0.get(), true);
    EXPECT_EQ(s.load(), 0);
  }

  {
    // now lets, see if we can interrupt the waiting thread by blocking
    std::promise<bool> p1;
    const auto task1 = [&s, &p1]() {
      try {
        s.consume_one({});
      } catch (const Interruptable::interrupt_set& e) {
        p1.set_exception(std::current_exception());
        return;
      } catch (const Interruptable::block_set& e) {
        p1.set_exception(std::current_exception());
        return;
      }
      p1.set_value_at_thread_exit(true);
    };

    std::thread t1(task1);
    auto f1 = p1.get_future();

    s.set_block(); // block-state set, next wake-up will cause a throw
    s.post_one(); //triggers a wake-up, before actual consumption, exception is raised
    t1.join();
    EXPECT_THROW(f1.get(), Interruptable::block_set);
    EXPECT_EQ(s.load(), 1);  //element was never consumed
  }


  {
    s.reset_block();

    bool end_exec{false};

    // now lets, see if we can interrupt the waiting thread
    std::promise<bool> p1;
    const auto task1 = [&s, &p1, &end_exec]() {
      try {
        s.consume_one([&]->bool{return end_exec;});
      } catch (const Interruptable::interrupt_set& e) {
        p1.set_exception(std::current_exception());
        return;
      } catch (const Interruptable::block_set& e) {
        p1.set_exception(std::current_exception());
        return;
      }
      p1.set_value_at_thread_exit(true);
    };

    std::thread t1(task1);
    std::this_thread::sleep_for(std::chrono::milliseconds(500));
    auto f1 = p1.get_future();

    end_exec = true;
    s.interrupt();
    t1.join();
    EXPECT_THROW(f1.get(), Interruptable::interrupt_set);
    EXPECT_EQ(s.load(), 0);
  }
}





TEST(Semaphore, interrupts) {
  Semaphore s;

  // interrupts that might be distributed among external (otherwise inaccessible)
  //  tasks
  auto i0 = s.get_interrupt();
  auto i1 = s.get_interrupt();

  EXPECT_FALSE(s.is_blocked());

  // set the trigger on i0
  EXPECT_FALSE(i0.is_triggered());
  i0.block();
  EXPECT_TRUE(i0.is_triggered());
  EXPECT_TRUE(s.is_blocked());
  i0.block();  // triggering an already triggered Interrupt is inert

  // reset the i0 trigger
  i0.unblock();
  EXPECT_FALSE(i0.is_triggered());
  EXPECT_FALSE(s.is_blocked());

  i0.unblock();  // resetting an already reset Interrupt is inert

  // trigger interrupt from multiple sites
  EXPECT_FALSE(s.is_blocked());
  i0.block();
  i1.block();
  EXPECT_TRUE(s.is_blocked());
  i0.unblock();
  EXPECT_TRUE(s.is_blocked());
  i1.unblock();
  EXPECT_FALSE(s.is_blocked());
}