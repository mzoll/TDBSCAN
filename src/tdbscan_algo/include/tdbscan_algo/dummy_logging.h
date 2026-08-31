//
// Created by netsu on 27/08/2026.
//

#ifndef TDBSCAN_DUMMY_LOGGING_H
#define TDBSCAN_DUMMY_LOGGING_H

#include <iostream>

void log_fatal(const std::string msg) {std::cout << "FATAL: " << msg << std::endl;};
void log_error(const std::string msg) {std::cout << "ERROR: " << msg << std::endl;};
void log_warn(const std::string msg) {std::cout << "WARN: " << msg << std::endl;};
void log_info(const std::string msg) {std::cout << "INFO: " << msg << std::endl;};
void log_debug(const std::string msg) {std::cout << "DEBUG: " << msg << std::endl;};
void log_trace(const std::string msg) {std::cout << "TRACE: " << msg << std::endl;};


#endif //TDBSCAN_DUMMY_LOGGING_H
