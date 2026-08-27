cmake_minimum_required (VERSION 3.11)

#include(cmake/StandardProjectSettings.cmake)

# User-configurable options
option(BUILD_EXAMPLES "Build example applications" ON)
option(BUILD_STATIC_LIB "Build static library only" OFF)
option(ENABLE_TIME_TRACE "Add -ftime-trace to Clang builds" OFF)

project (TDBSCAN
    VERSION 0.0.1
    LANGUAGES CXX
#    CXX_STANDARD 17  #current
#    CXX_STANDARD_REQUIRED ON
)
#
#set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
#set(CMAKE_COMPILE_WARNING_AS_ERROR ON)
#
## setup testing
#include(CTest)
#enable_testing()
#
#
#if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
#    set_property(GLOBAL PROPERTY USE_FOLDERS ON)
#endif()
#
#include(FetchContent)
#
#include_directories("${gtest_SOURCE_DIR}/include"
#        "${PROJECT_SOURCE_DIR}/include"
#)


add_subdirectory(src)

#add_subdirectory(test)