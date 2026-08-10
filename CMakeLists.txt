include(cmake/StandardProjectSettings.cmake)

cmake_minimum_required (VERSION 3.6)

project (TDBSCAN)

add_executable(TDBSCAN src/TDBSCAN.cpp
        lib/pure_algo/include/tdbscan_algo.h
        lib/pure_algo/include/tdbcluster.h
        lib/pure_algo/include/base_defs.h
        lib/pure_algo/include/absblib.h
        lib/pure_algo/include/tdbscan_algo.hh
        lib/pure_algo/include/connector.h
        lib/pure_algo/tdbcluster.cpp)
