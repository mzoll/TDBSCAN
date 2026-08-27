cmake_minimum_required (VERSION 3.6)

#include(cmake/StandardProjectSettings.cmake)

# User-configurable options
option(BUILD_EXAMPLES "Build example applications" ON)
option(BUILD_STATIC_LIB "Build static library only" OFF)
option(ENABLE_TIME_TRACE "Add -ftime-trace to Clang builds" OFF)

project (TDBSCAN
    VERSION 0.0.1
    LANGUAGES CXX
    CXX_STANDARD 14
    CXX_STANDARD_REQUIRED ON
)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

add_subdirectory(src)

add_subdirectory(test)