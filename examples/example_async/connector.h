//
// Created by netsu on 08/10/2026.
//

#ifndef TDBSCAN__EXASYNC__CONNECTOR_H
#define TDBSCAN__EXASYNC__CONNECTOR_H

#include "tdbscan_algo/tdbscan_algo.h"
#include "tdbscan_algo/common_blibs.h"


namespace exasync {

using namespace std;
using namespace tdbscan;

/**
 * Limits the Distance within which ScalarBlibs can causally connect
 */
class DistanceLimiter : public ConnectorSingle<ScalarBlib> {
protected:
  ScalarBlib::Ordinate_t::Distance_t maxDist_;

public:
  DistanceLimiter(const ScalarBlib::Ordinate_t::Distance_t maxDistance) : ConnectorSingle("DistanceLimiter"),
                                                                               maxDist_(maxDistance) {};

  [[nodiscard]] bool eval(const ScalarBlib &lhs, const ScalarBlib &rhs) const override {
    return fabs(lhs.distanceTo(rhs)) <= maxDist_;
  };
};


/**
 * Limits the Distance within which Blibs causally connect
 */
class TimeLimiter : public ConnectorSingle<ScalarBlib> {
private:
  ScalarBlib::Time_t::TimeDiff_t maxTimediff_;

public:
  TimeLimiter(const ScalarBlib::Time_t::TimeDiff_t maxTimeDiff)
  : ConnectorSingle("TimeLimiter"), maxTimediff_(maxTimeDiff) {};

  [[nodiscard]] bool eval(const ScalarBlib &lhs, const ScalarBlib &rhs) const final {
    return fabs(lhs.timeTo(rhs)) <= maxTimediff_;
  }
};


/// combine the Connectors into a ConnectorBlock
class LimitingConnector final : public ConnectorAssembly_AND<ScalarBlib> {
  const DistanceLimiter *const distance_limiter_;
  const TimeLimiter *const time_limiter_;

public:
  LimitingConnector(const double distance_lim, const double time_lim)
  : distance_limiter_(
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

} // namespace exasync

#endif //TDBSCAN__EXASYNC__CONNECTOR_H
