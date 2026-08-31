//
// Created by netsu on 30/08/2026.
//


#include <gtest/gtest.h>
#include "tdbscan_algo/connector.h"

#include "tdbscan_algo/common_defs.h"

using namespace tdbscan;


// having implemented all this stuff, lets whip up a connector, which is just an distance connector
class DistanceLimiter final : public tdbscan::ConnectorSingle<Blib4d> {
public:
  Blib4d::Ordinate_t::Distance_t maxDist_;
  DistanceLimiter(const Blib4d::Ordinate_t::Distance_t maxDistance) : ConnectorSingle("DistConnector"), maxDist_(maxDistance) {};

  bool eval(const Blib4d& lhs, const Blib4d& rhs) const {return lhs.getDistance(rhs) <= maxDist_;};
};

// make one connector which just connects to max time-diff
class TimeLimiter final : public ConnectorSingle<Blib4d> {
public:
  Blib4d::Time_t::TimeDiff_t maxTimediff_;
  explicit TimeLimiter(const Blib4d::Time_t::TimeDiff_t maxTimeDiff) : ConnectorSingle("DistConnector"), maxTimediff_(maxTimeDiff) {};

  bool eval(const Blib4d& lhs, const Blib4d& rhs) const {return rhs.timeDiff(lhs) <= maxTimediff_;};
};

// Demonstrate some basic assertions.
TEST(ConnectorTest, DummyConnector) {
  auto distLimiter_ = new DistanceLimiter(20.);

  //evaluate
  distLimiter_->eval(Blib4d({0,0,0},0), Blib4d({0,0,0},0));
  delete distLimiter_;
}


class LimitingConnector final : public ConnectorBlock<Blib4d> {};


TEST(ConnectorTest, DummyConnectorBlock) {

  auto distLimiter_ = new DistanceLimiter(20.);
  auto timeLimiter_ = new TimeLimiter(20.);
  auto limcon = new LimitingConnector();
  limcon->addConnector(distLimiter_);
  limcon->addConnector(timeLimiter_);

  limcon->eval(Blib4d({0,0,0},0), Blib4d({0,0,0},0));
};



