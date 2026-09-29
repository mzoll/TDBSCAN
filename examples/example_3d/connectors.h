//
// Created by netsu on 27/09/2026.
//

#ifndef TDBSCAN__EXAMPLE_3D__CONNECTORS_H
#define TDBSCAN__EXAMPLE_3D__CONNECTORS_H

#include "blib.h"
#include <tdbscan_algo/connector.h>


namespace ex3d {
using namespace tdbscan;

// define some Limiters
class DistanceLimiter_ : public ConnectorSingle<Blib3d> {
public:
  Blib3d::Ordinate_t::Distance_t maxDist_;

  DistanceLimiter_(const Blib3d::Ordinate_t::Distance_t maxDistance) : ConnectorSingle("DistConnector"),
                                                                       maxDist_(maxDistance) {};

  bool eval(const Blib3d &lhs, const Blib3d &rhs) const { return lhs.distanceTo(rhs) <= maxDist_; };
};


/// Extends the DistanceLimiter_ to SBlibWithTrace
class DistanceLimiter : public DistanceLimiter_, public ConnectorSingle<Blib3dWithTrace> {
public:
  DistanceLimiter(const Blib3dWithTrace::Ordinate_t::Distance_t maxDistance)
    : DistanceLimiter_(maxDistance),
      ConnectorSingle<Blib3dWithTrace>(DistanceLimiter_::name_) {};

  [[nodiscard]] inline
  bool eval(const Blib3dWithTrace &lhs, const Blib3dWithTrace &rhs) const {
    return DistanceLimiter_::eval(lhs, rhs);
  };
};

// make one connector which just connects to max time-diff
class TimeLimiter_ : public ConnectorSingle<Blib3d> {
public:
  Blib3d::Time_t::TimeDiff_t maxTimediff_;

  explicit TimeLimiter_(const Blib3d::Time_t::TimeDiff_t maxTimeDiff) : ConnectorSingle("DistConnector"),
                                                                        maxTimediff_(maxTimeDiff) {};

  bool eval(const Blib3d &lhs, const Blib3d &rhs) const { return rhs.timeTo(lhs) <= maxTimediff_; };
};


/// Extends the DistanceLimiter_ to SBlibWithTrace
class TimeLimiter : public TimeLimiter_, public ConnectorSingle<Blib3dWithTrace> {
public:
  TimeLimiter(const Blib3dWithTrace::Ordinate_t::Distance_t maxDistance)
    : TimeLimiter_(maxDistance),
      ConnectorSingle<Blib3dWithTrace>(TimeLimiter_::name_) {};

  [[nodiscard]] inline
  bool eval(const Blib3dWithTrace &lhs, const Blib3dWithTrace &rhs) const {
    return TimeLimiter_::eval(lhs, rhs);
  };
};


// combine the Connectors into a ConnectorBlock
class LimitingConnector final : public ConnectorAssembly_AND<Blib3dWithTrace> {
  const DistanceLimiter *const distance_limiter_;
  const TimeLimiter *const time_limiter_;

public:
  LimitingConnector(const double distance_lim, const double time_lim) : distance_limiter_(
                                                                          new DistanceLimiter(distance_lim)),
                                                                        time_limiter_(new TimeLimiter(time_lim)) {
    addConnector(distance_limiter_);
    addConnector(time_limiter_);
  };

  ~LimitingConnector() override {
    delete distance_limiter_;
    delete time_limiter_;
  };
};
} // namespace ex3d

#endif //TDBSCAN__EXAMPLE_3D__CONNECTORS_H
