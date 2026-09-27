//
// Created by netsu on 27/09/2026.
//

#ifndef TDBSCAN__EXAMPLE_3D__BLIB_H
#define TDBSCAN__EXAMPLE_3D__BLIB_H


#include "tdbscan_algo/common_defs.h"

/**
 * A twist to the Blib3d class that allows to keep track of its source, either Noise or Signal
 */
class Blib3dWithTrace final : public Blib3d {
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
  Blib3dWithTrace& mark(const Origin o) { origin_ = o; return *this; };


  /**
   * Fully qualified constructor
   * @param ord the Ordinate
   * @param t time of the blib
   * @param o the origin of the Blib, aka SIGNAL, or NOISE
   */
  Blib3dWithTrace(
    const Blib3dWithTrace::Ordinate_t& ord,
    const Blib3dWithTrace::Time_t& t,
    const Origin o ) :
  Blib3d(ord, t) ,origin_(o) {};

  inline
  bool operator==(const Blib3dWithTrace& rhs) const {
    return this->origin_ == rhs.origin_ && static_cast<Blib3d>(*this) == static_cast<Blib3d>(rhs);
  };

public:
  friend
  std::ostream& operator<< ( std::ostream& outs, const Blib3dWithTrace & st);
};


/// custom ostream output; nicely formats the class as a string
inline
std::ostream& operator<< ( std::ostream& os, const Blib3dWithTrace & sblib_trace) {
  return os <<
    Blib3dWithTrace::origin_tostr(sblib_trace.origin_) << "::" <<
      static_cast<Blib3d>(sblib_trace);
};

/// custom formatter for Blib3dWithTrace. Used in conjunction with `std::format` or `fmt::format`
template <>
struct std::formatter<Blib3dWithTrace> {
  constexpr auto parse(std::format_parse_context & ctx) {return ctx.begin();}

  auto format(const Blib3dWithTrace& sb, std::format_context &ctx) const {
    return std::format_to(ctx.out(), "[ord:{}, time:{}]::{}",
      sb.getOrdinate(),
      static_cast<double>(sb.getTime()),
      Blib3dWithTrace::origin_tostr(sb.origin_));
  };
};


#endif //TDBSCAN__EXAMPLE_3D__BLIB_H
