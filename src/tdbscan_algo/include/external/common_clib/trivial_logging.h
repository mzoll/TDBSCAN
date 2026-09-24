//
// Created by mzoll on 02.03.21.
//
// code example from here
// https://www.scalyr.com/blog/getting-started-quickly-c++-logging also here a
// customization
// https://stackoverflow.com/questions/53744798/boost-boost-log-sev-with-custom-attributes

#ifndef COMMONCLIB_TRIVIAL_LOGGING_H
#define COMMONCLIB_TRIVIAL_LOGGING_H

#include <boost/log/core.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/utility/setup/common_attributes.hpp>
#include <boost/log/utility/setup/file.hpp>

// if set true DEBUG and TRACE logs will be optimized out
#define OPTIMIZE_RELEASE_COMPILE 1

#if (OPTIMIZED_RELEASE_COMPILE + NDEBUG >= 2)
#define RELEASE_OPT
#endif

// System Log macros.
// TRACE < DEBUG < INFO < WARN < ERROR < FATAL
// DEBUG (with optimzations)

// --- make a convenience formatter with fmt
#include <fmt/format.h>
#include <fmt/ostream.h>

#ifdef RELEASE_OPT
#define LOG_TRACE(...) ;
#define LOG_DEBUG(...) ;
#else
#define LOG_TRACE(...) \
  BOOST_LOG_TRIVIAL(trace) << fmt::format(__VA_ARGS__);
#define LOG_DEBUG(...) \
  BOOST_LOG_TRIVIAL(debug) << fmt::format(__VA_ARGS__);
#endif

#define LOG_INFO(...) \
  BOOST_LOG_TRIVIAL(info) << fmt::format(__VA_ARGS__);
#define LOG_WARNING(...) \
  BOOST_LOG_TRIVIAL(warning) << fmt::format(__VA_ARGS__);
#define LOG_ERROR(...) \
  BOOST_LOG_TRIVIAL(error) << fmt::format(__VA_ARGS__);
#define LOG_FATAL(...) \
  BOOST_LOG_TRIVIAL(fatal) << fmt::format(__VA_ARGS__);

#endif  // COMMONCLIB_TRIVIAL_LOGGING_H