add_library(tdbscan_algo)

target_include_directories(tdbscan_algo
        PUBLIC include
        PRIVATE src
)

target_link_libraries(tdbscan_algo

)

target_sources(hello_world
        PRIVATE
        src/TDBSCAN.cpp)

target_compile_features(hello_world PUBLIC cxx_std_17)

set_target_properties(hello_world PROPERTIES
        CXX_STANDARD_REQUIRED ON
        CXX_EXTENSIONS OFF
)

# /// IDEA ///
# https://dev.to/pigeoncodeur/cmake-for-complex-projects-part-1-building-a-c-game-engine-from-scratch-for-desktop-and-ied
# Modern CMake magic - precompiled headers for faster builds
# target_precompile_headers(ColumbaEngine PUBLIC src/Engine/stdafx.h)

target_sources(pure_algo
        PUBLIC
        src/pure_algo/include/tdbscan_algo.h
        src/pure_algo/include/tdbcluster.h
        src/pure_algo/include/base_defs.h
        src/pure_algo/include/absblib.h
        src/pure_algo/include/tdbscan_algo.hh
        src/pure_algo/include/connector.h
        PRIVATE
        src/pure_algo/tdbcluster.cpp
)


target_sources(Testsuite
    PRIVATE
    test/pure_algo/test_absblib.cpp
)
