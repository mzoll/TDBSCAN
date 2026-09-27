//
// Created by netsu on 27/09/2026.
//

#ifndef TDBSCAN__EXAMPLE_3D__CONNECTORS_H
#define TDBSCAN__EXAMPLE_3D__CONNECTORS_H

#include "blib.h"
#include <tdbscan_algo/connector.h>


// define some Limiters
class DistanceLimiter final : public tdbscan::ConnectorSingle<Blib3d> {
public:
  Blib3d::Ordinate_t::Distance_t maxDist_;
  DistanceLimiter(const Blib3d::Ordinate_t::Distance_t maxDistance) : ConnectorSingle("DistConnector"), maxDist_(maxDistance) {};

  bool eval(const Blib3d& lhs, const Blib3d& rhs) const {return lhs.distanceTo(rhs) <= maxDist_;};
};

// make one connector which just connects to max time-diff
class TimeLimiter final : public tdbscan::ConnectorSingle<Blib3d> {
public:
  Blib3d::Time_t::TimeDiff_t maxTimediff_;
  explicit TimeLimiter(const Blib3d::Time_t::TimeDiff_t maxTimeDiff) : ConnectorSingle("DistConnector"), maxTimediff_(maxTimeDiff) {};

  bool eval(const Blib3d& lhs, const Blib3d& rhs) const {return rhs.timeTo(lhs) <= maxTimediff_;};
};

// combine the Connectors into a ConnectorBlock
class LimitingConnector final : public tdbscan::ConnectorBlock<Blib3d> {
  const DistanceLimiter* const distance_limiter_;
  const TimeLimiter* const time_limiter_;
public:
  LimitingConnector(const double distance_lim, const double time_lim) :
    distance_limiter_(new DistanceLimiter(distance_lim)),
    time_limiter_(new TimeLimiter(time_lim)) {
    addConnector(distance_limiter_);
    addConnector(time_limiter_);
  };

  ~LimitingConnector() override {
    delete distance_limiter_;
    delete time_limiter_;
  };
};

#endif //TDBSCAN__EXAMPLE_3D__CONNECTORS_H
