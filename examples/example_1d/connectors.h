//
// Created by netsu on 27/09/2026.
//

#ifndef TDBSCAN__EX1D__CONNECTORS_H
#define TDBSCAN__EX1D__CONNECTORS_H


#include "tdbscan_algo/connector.h"
#include "blib.h"

namespace ex1d {
using namespace tdbscan;
// ===================== Define some Limiters ====================
// We create two simple limiters: One to limit distance and one to Limit time between Blibs
// Both of them are Connectors and compare one Blib to another.
// The only thing we need to implement is the the `bool eval(Blib rhs, Blib rhs)´ function.
// In the end we combine both Limiters into a single `ConnectorBlock`, which than can be used in the TDBscan-algorithm
// ========================================================

/**
 * Limits the Distance within which ScalarBlibs can causally connect
 */
class DistanceLimiter_ : public ConnectorSingle<ScalarBlib> {
protected:
  SBlibWithTrace::Ordinate_t::Distance_t maxDist_;

public:
  DistanceLimiter_(const SBlibWithTrace::Ordinate_t::Distance_t maxDistance) : ConnectorSingle("DistanceLimiter"),
                                                                               maxDist_(maxDistance) {};

  [[nodiscard]] bool eval(const ScalarBlib &lhs, const ScalarBlib &rhs) const override {
    return fabs(lhs.distanceTo(rhs)) <= maxDist_;
  };
};

/// Extends the DistanceLimiter_ to SBlibWithTrace
class DistanceLimiter : public DistanceLimiter_, public ConnectorSingle<SBlibWithTrace> {
public:
  DistanceLimiter(const SBlibWithTrace::Ordinate_t::Distance_t maxDistance)
    : DistanceLimiter_(maxDistance),
      ConnectorSingle<SBlibWithTrace>(DistanceLimiter_::name_) {};

  [[nodiscard]] inline
  bool eval(const SBlibWithTrace &lhs, const SBlibWithTrace &rhs) const {
    return DistanceLimiter_::eval(lhs, rhs);
  };
};

/**
 * Limits the Distance within which Blibs causally connect
 */
class TimeLimiter_ : public ConnectorSingle<ScalarBlib> {
private:
  SBlibWithTrace::Time_t::TimeDiff_t maxTimediff_;

public:
  TimeLimiter_(const SBlibWithTrace::Time_t::TimeDiff_t maxTimeDiff) : ConnectorSingle("TimeLimiter"),
                                                                       maxTimediff_(maxTimeDiff) {};

  [[nodiscard]] bool eval(const ScalarBlib &lhs, const ScalarBlib &rhs) const final {
    return fabs(lhs.timeTo(rhs)) <= maxTimediff_;
  }
};

/// Extends the TimeLimiter_ to SBlibWithTrace
class TimeLimiter : public TimeLimiter_, public ConnectorSingle<SBlibWithTrace> {
public:
  TimeLimiter(const SBlibWithTrace::Ordinate_t::Distance_t maxDistance)
    : TimeLimiter_(maxDistance),
      ConnectorSingle<SBlibWithTrace>(TimeLimiter_::name_) {};

  [[nodiscard]] inline
  bool eval(const SBlibWithTrace &lhs, const SBlibWithTrace &rhs) const {
    return TimeLimiter_::eval(lhs, rhs);
  };
};


/**
 * Assumes that hits are intrinsically caused by a source that moves with a certain inertia.
 */
class InertiaConnector_ : public ConnectorSingle<ScalarBlib> {
private:
  const double inertia_;
  const double tollerance_abs_;

public:
  explicit InertiaConnector_(const double inertia, const double tollerance_abs) : ConnectorSingle("InertiaConnector"),
    inertia_(inertia), tollerance_abs_(tollerance_abs) {};

  [[nodiscard]] bool eval(const ScalarBlib &lhs, const ScalarBlib &rhs) const final {
    return fabs(lhs.timeTo(rhs) * inertia_ - lhs.distanceTo(rhs)) <= tollerance_abs_;
  };
};


class InertiaConnector : public InertiaConnector_, public ConnectorSingle<SBlibWithTrace> {
  InertiaConnector(const double inertia, const double tollerance_abs) : InertiaConnector_(inertia, tollerance_abs),
                                                                        ConnectorSingle<SBlibWithTrace>(
                                                                          ConnectorSingle<ScalarBlib>::name_) {};
};


/// combine the Connectors into a ConnectorBlock
class LimitingConnector final : public ConnectorAssembly_AND<SBlibWithTrace> {
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
} // namespace ex1d


#endif //TDBSCAN__EX1D__CONNECTORS_H
