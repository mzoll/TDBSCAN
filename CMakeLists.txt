include(cmake/StandardProjectSettings.cmake)

cmake_minimum_required (VERSION 3.6)

project (TDBSCAN
    VERSION 0.0.1
    LANGUAGES CXX
    CXX_STANDARD 14
    CXX_STANDARD_REQUIRED ON
)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

add_subdirectory(src)

add_subdirectory(test)