//
// Created by mzoll on 31/08/2026.
//


#include <condition_variable>
#include <thread>
#include <iostream>
#include <atomic>
#include <mutex>
#include <functional>

#include <chrono>
using namespace std::chrono_literals;

using namespace std;


struct ConsumptionEnded : std::runtime_error {
  ConsumptionEnded() : std::runtime_error("ConsumptionEnded") {};
};

class Sem {
public:
  condition_variable cond;
private:
  std::atomic<int> counter{4};
  std::mutex mtx;
public:
  bool end_consumption{false};
public:
  int consume(std::function<void()> fkt, std::chrono::duration<int> timeout = std::chrono::seconds(1)) {
    cout << "Consume requested: " << counter.load() << endl;
    if (end_consumption)
      throw std::runtime_error("consumption is interrupted");
    std::unique_lock<mutex> lock(mtx);
    if (counter.load() == 0) {
      //cond.wait(lock,[&] {
      cond.wait_for(lock, timeout, [&] {
        cout << "inner predicate" <<endl;
        if (end_consumption)
          throw ConsumptionEnded();
        fkt();
        return counter > 0;
      });
      cout << "Thinking there is something to consume" << endl;
    }

  cout << "Ready to consume; counter : " << counter.load() << endl;;
  --counter;
  if (counter.load() < 0)
    throw std::logic_error("Counter is negative !!! ");
  return 1;
  }
};


struct Disruption : std::runtime_error {
  Disruption() : std::runtime_error("Disruption") {};
};


int main() {
  Sem sem;

  bool is_disrupted{false};

  // auto disruption = [&]() {
  //   cout << "probing disruption" << endl;
  //   if (is_disrupted) {
  //     cout << "throwing disruption" << endl;
  //     throw Disruption();
  //   }
  // };

  auto disruption = []() {
    static int count_down = 4;
    cout << "probing disruption" << endl;
    if (count_down <= 0) {
      cout << "throwing disruption" << endl;
      throw Disruption();
    }
    --count_down;
  };

  thread t0([&]() {
    while (true) {
      try{
        cout << ">> Consume One" <<endl;
        sem.consume(disruption, std::chrono::seconds(1));
        cout << "<< Consumed One" <<endl;
      }
      catch (ConsumptionEnded& e) {
        cout << "resurfacing: internal" <<endl;
        break;
      }
      catch (Disruption& e) {
        cout << "resurfacing: disruption" <<endl;
        break;
      }
    }
  });

  std::this_thread::sleep_for(3s);
  //cout << " xxx trigger" <<endl;
  //is_disrupted = true;
  //sem.is_blocked = true;
  //sem.cond.notify_all();

  t0.join();

  return 0;
}
