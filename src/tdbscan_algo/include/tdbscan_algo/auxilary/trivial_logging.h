//
// Created by mzoll on 02.03.21.
//
// code example from here
// https://www.scalyr.com/blog/getting-started-quickly-c++-logging also here a
// customization
// https://stackoverflow.com/questions/53744798/boost-boost-log-sev-with-custom-attributes

#ifndef COMMONCLIB_TRIVIAL_LOGGING_H
#define COMMONCLIB_TRIVIAL_LOGGING_H

#include <iostream>
#include <format>

// if set true DEBUG and TRACE logs will be optimized out
#define OPTIMIZE_RELEASE_COMPILE 1

#if (OPTIMIZED_RELEASE_COMPILE + NDEBUG >= 2)
#define RELEASE_OPT
#endif

// System Log macros.
// TRACE < DEBUG < INFO < WARN < ERROR < FATAL
// DEBUG (with optimisations)

#define RELEASE_OPT

#ifdef RELEASE_OPT
#define LOG_TRACE(...) ;
#define LOG_DEBUG(...) ;
#else
#define LOG_TRACE(...) \
  std::cout << "TRACE: " << std::format(__VA_ARGS__) << std::endl;
#define LOG_DEBUG(...) \
  std::cout << "DEBUG: " << std::format(__VA_ARGS__) << std::endl;
#endif

#define LOG_INFO(...) \
  std::cout << "INFO: " << std::format(__VA_ARGS__) << std::endl;
#define LOG_WARNING(...) \
  std::cout << "WARNING: " << std::format(__VA_ARGS__) << std::endl;
#define LOG_ERROR(...) \
  std::cout << "ERROR: " << std::format(__VA_ARGS__) << std::endl;
#define LOG_FATAL(...) \
  std::cout << "FATAL: " << std::format(__VA_ARGS__) << std::endl;

#endif  // COMMONCLIB_TRIVIAL_LOGGING_H