//
// Created by marcel on 23.02.21.
//


#include <gtest/gtest.h>

#include "external/common_clib/InterruptableQueue.hpp"

#include <chrono>
#include <exception>
#include <future>
#include <thread>


using namespace std;
using namespace common_clib;
using namespace common_clib::threading;


TEST(InterruptableQueue, serial_access) {


  InterruptableQueue<int> queue;

  EXPECT_FALSE(queue.is_closed());

  EXPECT_EQ(queue.size(), 0);
  EXPECT_TRUE(queue.empty());

  int e_counter = 0;
  const int element_0 = ++e_counter;
  const int element_1 = ++e_counter;

  {
    // single element
    EXPECT_TRUE(queue.empty());
    EXPECT_NO_THROW(queue.push(element_0));
    EXPECT_NO_THROW(queue.push(element_1));
    EXPECT_EQ(queue.size(), 2);
    EXPECT_FALSE(queue.empty());
    int retrieval;
    EXPECT_NO_THROW(retrieval = queue.pop());
    EXPECT_EQ(queue.size(), 1);
    EXPECT_EQ(retrieval, element_0);  // elements arrive in the right order: FIFO
    EXPECT_NO_THROW(retrieval = queue.pop());
    EXPECT_EQ(queue.size(), 0);
    EXPECT_EQ(retrieval, element_1);  // elements arrive in the right order: FIFO
  }

  {
    //many elements
    EXPECT_NO_THROW(queue.push(element_0));
    EXPECT_NO_THROW(queue.push(element_1));
    EXPECT_EQ(queue.size(), 2);
    std::vector<int> retrieval;
    EXPECT_NO_THROW(retrieval = queue.exhaust());
    EXPECT_EQ(queue.size(), 0);
    EXPECT_EQ(retrieval[0], element_0);
    EXPECT_EQ(retrieval[1], element_1);
  }
}

TEST(InterruptableQueue, serial_access_blocking) {

  int _;

  InterruptableQueue<int> queue;

  EXPECT_FALSE(queue.is_closed());

  EXPECT_NO_THROW(queue.push(42));
  EXPECT_EQ(queue.size(), 1);

  // block the outlet,
  EXPECT_NO_THROW(queue.close_outlet());
  EXPECT_TRUE(queue.is_closed());
  EXPECT_NO_THROW(queue.push(42));
  EXPECT_EQ(queue.size(), 2);
  EXPECT_ANY_THROW(_= queue.pop());

  // block the inlet too
  EXPECT_NO_THROW(queue.close_inlet());
  EXPECT_TRUE(queue.is_closed());
  EXPECT_ANY_THROW(queue.push(42));
  EXPECT_EQ(queue.size(), 2);
  EXPECT_ANY_THROW(_= queue.pop());

  //release the outlet
  EXPECT_NO_THROW(queue.reopen_outlet());
  EXPECT_TRUE(queue.is_closed());
  EXPECT_ANY_THROW(queue.push(42));
  EXPECT_EQ(queue.size(), 2);
  EXPECT_NO_THROW(_= queue.pop());
  EXPECT_EQ(queue.size(), 1);

  //release the inlet
  EXPECT_NO_THROW(queue.reopen_inlet());
  EXPECT_FALSE(queue.is_closed());
  EXPECT_NO_THROW(queue.push(42));
  EXPECT_EQ(queue.size(), 2);
  EXPECT_NO_THROW(_= queue.pop());
  EXPECT_EQ(queue.size(), 1);
}

TEST(InterruptableQueue, parallel_access) {
  InterruptableQueue<int> q;
  {
    std::promise<int> p0;

    const auto task0 = [&q, &p0]() {
      int retrieval;
      try {
        retrieval = q.pop();
      } catch (const InterruptableQueue<int>::InterruptIsSet& e) {
        p0.set_exception(std::current_exception());
        return;
      }
      p0.set_value_at_thread_exit(retrieval);
    };


    std::thread t0(task0);
    std::this_thread::sleep_for(std::chrono::milliseconds(500));
    auto f0 = p0.get_future();
    f0.wait_for(std::chrono::milliseconds(100));

    q.push(42);  //push something to the queue
    t0.join();  //eventually t0 will be notified
    EXPECT_EQ(f0.get(), 42);
    EXPECT_EQ(q.size(), 0);
  }

  {
    // now lets, see if we can set an interrupt on the waiting
    std::promise<bool> p1;
    const auto task1 = [&q, &p1]() {
      try {
        p1.set_value_at_thread_exit( q.pop() );
      } catch (const InterruptableQueue<int>::InterruptIsSet& e) {
        p1.set_exception(std::current_exception());
        return;
      }
    };

    std::thread t1(task1);
    std::this_thread::sleep_for(std::chrono::milliseconds(500));
    auto f1 = p1.get_future();

    q.close_outlet();  // the outlet is now blocked
    q.push(42);  // push a value; if everything works correctly, it is never retrieved
    t1.join();  // by blocking the outlet, the waiting of t1 was interrupted
    EXPECT_THROW(f1.get(), InterruptableQueue<int>::InterruptIsSet);
    EXPECT_TRUE(q.size() == 1);
  }
}

TEST(InterruptableQueue, interrupts) {
  //FIXME write this test using the same mechanism as in test_semaphore
  // same as in
}