//
// Created by mzoll on 30/08/2026.
//

#ifndef TDBSCAN__COMMON_BLIBS_HH
#define TDBSCAN__COMMON_BLIBS_HH

#include <format>

#include "common_defs.h"

inline
std::ostream&
operator<<(std::ostream& os, const tdbscan::Blib3d & b4d) {
  return os << std::format("Blib3d(ord:{}, time:{})", b4d.getOrdinate(), b4d.getTime());
};

template <>
struct std::formatter<tdbscan::Blib3d> {
  constexpr auto parse(std::format_parse_context & ctx) {return ctx.begin();}

  auto format(const tdbscan::Blib3d& sb, std::format_context &ctx) const {
    return std::format_to(ctx.out(), "[ord:{}, time:{}]", sb.getOrdinate(), static_cast<double>(sb.getTime()));
  };
};

inline
std::ostream&
operator<<(std::ostream& os, const tdbscan::ScalarBlib & sb) {
  return os << std::format("[ord: {}, time: {}]", double(sb.getOrdinate()), double(sb.getTime()));
};

template <>
struct std::formatter<tdbscan::ScalarBlib> {
  constexpr auto parse(std::format_parse_context & ctx) {return ctx.begin();}

  auto format(const tdbscan::ScalarBlib& sb, std::format_context &ctx) const {
    return std::format_to(ctx.out(), "[ord:{}, time:{}]", static_cast<double>(sb.getOrdinate()), static_cast<double>(sb.getTime()) );
  };
};

#endif //TDBSCAN__COMMON_BLIBS_HH
