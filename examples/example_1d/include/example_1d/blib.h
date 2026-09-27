//
// Created by netsu on 27/09/2026.
//

#ifndef TDBSCAN__EX1D__BLIB_H
#define TDBSCAN__EX1D__BLIB_H

#include "tdbscan_algo/common_defs.h"

/**
 * A twist to the ScalarBlib class that allows to keep track of its source, either Noise or Signal
 */
class SBlibWithTrace final : public ScalarBlib {
public:
  /// denotes the origin of the Blib; either Noise or Signal
  enum Origin {
    UNKNOWN = 0,
    NOISE = 20,
    SIGNAL = 100,
  } origin_{UNKNOWN};

  ///convenience to convert `Origin` into a string for human reading
  inline
  static std::string origin_tostr(const Origin o) {
    switch (o) {
      case UNKNOWN: return "UNKNOWN";
      case NOISE: return "NOISE";
      case SIGNAL: return "SIGNAL";
      default: throw std::invalid_argument("Value_error");
    }

  }

  ///mark this Blib as to stem either from a Noise or Signal source
  SBlibWithTrace& mark(const Origin o) { origin_ = o; return *this; };


  /**
   * Fully qualified constructor
   * @param ord the Ordinate
   * @param t time of the blib
   * @param o the origin of the Blib, aka SIGNAL, or NOISE
   */
  SBlibWithTrace(
    const SBlibWithTrace::Ordinate_t& ord,
    const SBlibWithTrace::Time_t& t,
    const Origin o ) :
  ScalarBlib(ord, t) ,origin_(o) {};

  inline
  bool operator==(const SBlibWithTrace& rhs) const {
    return this->origin_ == rhs.origin_ && static_cast<ScalarBlib>(*this) == static_cast<ScalarBlib>(rhs);
  };

public:
  friend
  std::ostream& operator<< ( std::ostream& outs, const SBlibWithTrace & st);
};


/// custom ostream output; nicely formats the class as a string
inline
std::ostream& operator<< ( std::ostream& os, const SBlibWithTrace & sblib_trace) {
  return os <<
    SBlibWithTrace::origin_tostr(sblib_trace.origin_) << "::" <<
      static_cast<ScalarBlib>(sblib_trace);
};

/// custom formatter for SBlibWithTrace. Used in conjunction with `std::format` or `fmt::format`
template <>
struct std::formatter<SBlibWithTrace> {
  constexpr auto parse(std::format_parse_context & ctx) {return ctx.begin();}

  auto format(const SBlibWithTrace& sb, std::format_context &ctx) const {
    return std::format_to(ctx.out(), "[ord:{}, time:{}]::{}",
      static_cast<double>(sb.getOrdinate()),
      static_cast<double>(sb.getTime()),
      SBlibWithTrace::origin_tostr(sb.origin_));
  };
};

#endif //TDBSCAN__EX1D__BLIB_H
