//
// Created by mzoll on 02.03.21.
//
// code example from here
// https://www.scalyr.com/blog/getting-started-quickly-c++-logging also here a
// customization
// https://stackoverflow.com/questions/53744798/boost-boost-log-sev-with-custom-attributes

#ifndef COMMONCLIB_TRIVIAL_LOGGING_H
#define COMMONCLIB_TRIVIAL_LOGGING_H

#include <format>
#include <exception>

// if set true DEBUG and TRACE logs will be optimized out
#define OPTIMIZE_RELEASE_COMPILE 1

#if (OPTIMIZED_RELEASE_COMPILE + NDEBUG >= 2)
#define RELEASE_OPT
#endif

#ifdef USE_BOOST_LOGGING
  // use the boost trivial logging
  #include <boost/log/core.hpp>
  #include <boost/log/expressions.hpp>
  #include <boost/log/trivial.hpp>
  #include <boost/log/utility/setup/common_attributes.hpp>
  #include <boost/log/utility/setup/file.hpp>

  #define LOG_FKT(log_level, ...) \
    BOOST_LOG_TRIVIAL(log_level) << std::format(__VA_ARGS__);
#else
  #include <iostream>

  inline
  std::string runtime_fallback(const std::string& x) { return x; };
  constexpr std::string _to_logstr(const std::string& llevel) {
    if (llevel == "trace")
        return "TRACE";
    if (llevel == "debug")
        return "DEBUG";
    if (llevel == "info")
        return "INFO";
    if (llevel == "warning")
        return "warning";
    if (llevel == "error")
        return "ERROR";
    if (llevel == "fatal")
        return "FATAL";
    throw std::domain_error(llevel);
  };
  #define LOG_FKT(log_level, ...) \
    std::cout << _to_logstr(#log_level) << ": " << std::format(__VA_ARGS__) << std::endl;
#endif

#define RELEASE_OPT

// System Log macros.
// TRACE < DEBUG < INFO < WARN < ERROR < FATAL
// DEBUG (with optimisations)

#ifdef RELEASE_OPT
#define LOG_TRACE(...) ;
#define LOG_DEBUG(...) ;
#else
#define LOG_TRACE(...) \
  LOG_FKT(trace, __VA_ARGS__)
#define LOG_DEBUG(...) \
  LOG_FKT(debug, __VA_ARGS__)
#endif

#define LOG_INFO(...) \
  LOG_FKT(info, __VA_ARGS__)
#define LOG_WARNING(...) \
  LOG_FKT(warning, __VA_ARGS__)
#define LOG_ERROR(...) \
  LOG_FKT(error, __VA_ARGS__)
#define LOG_FATAL(...) \
  LOG_FKT(fatal, __VA_ARGS__)

#endif  // COMMONCLIB_TRIVIAL_LOGGING_H