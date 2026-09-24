//
// Created by netsu on 27/08/2026.
//

#ifndef TDBSCAN_DUMMY_LOGGING_H
#define TDBSCAN_DUMMY_LOGGING_H

#include <iostream>
#include <format>
#include <sstream>

inline void log_fatal(const std::string msg) {std::cout << "FATAL: " << msg << std::endl;};
inline void log_error(const std::string msg) {std::cout << "ERROR: " << msg << std::endl;};
inline void log_warn(const std::string msg) {std::cout << "WARN: " << msg << std::endl;};
inline void log_info(const std::string msg) {std::cout << "INFO: " << msg << std::endl;};
inline void log_debug(const std::string msg) {std::cout << "DEBUG: " << msg << std::endl;};
inline void log_trace(const std::string msg) {std::cout << "TRACE: " << msg << std::endl;};

inline void log_fatal(const std::ostringstream& msg) {std::cout << "FATAL: " << msg.str() << std::endl;};
inline void log_error(const std::ostringstream& msg) {std::cout << "ERROR: " << msg.str() << std::endl;};
inline void log_warn(const std::ostringstream& msg) {std::cout << "WARN: " << msg.str() << std::endl;};
inline void log_info(const std::ostringstream& msg) {std::cout << "INFO: " << msg.str() << std::endl;};
inline void log_debug(const std::ostringstream& msg) {std::cout << "DEBUG: " << msg.str() << std::endl;};
inline void log_trace(const std::ostringstream& msg) {std::cout << "TRACE: " << msg.str() << std::endl;};

//
// template <class... Args>
// void log_fatal(const std::string& fmt_string, Args &&... params)
// {
//   const std::string msg = std::format(fmt_string, std::forward<Args>(params)...);
//   std::cout << "FATAL: " << msg << '\n';
// }
//
// template <class... Args>
// void log_error(const std::string& fmt_string, Args &&... params)
// {
//   const std::string msg = std::format(fmt_string, std::forward<Args>(params)...);
//   std::cout << "ERROR: " << msg << '\n';
// }
//
// template <class... Args>
// void log_warn(const std::string& fmt_string, Args &&... params)
// {
//   const std::string msg = std::format(fmt_string, std::forward<Args>(params)...);
//   std::cout << "WARN: " << msg << '\n';
// }
//
// template <class... Args>
// void log_info(const std::string& fmt_string, Args &&... params)
// {
//   const std::string msg = std::format(fmt_string, std::forward<Args>(params)...);
//   std::cout << "INFO: " << msg << '\n';
// }
//
// template <class... Args>
// void log_debug(const std::string& fmt_string, Args &&... params)
// {
//   const std::string msg = std::format(fmt_string, std::forward<Args>(params)...);
//   std::cout << "DEBUG: " << msg << '\n';
// }
//
// template <class... Args>
// void log_trace(const std::string& fmt_string, Args &&... params)
// {
//   const std::string msg = std::format(fmt_string, std::forward<Args>(params)...);
//   std::cout << "TRACE: " << msg << '\n';
// }


#endif //TDBSCAN_DUMMY_LOGGING_H
