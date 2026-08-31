add_library(tdbscan_algo
        src/include_all.cpp
        )

target_include_directories(tdbscan_algo
        PUBLIC include
        PRIVATE src
)


# /// IDEA ///
# https://dev.to/pigeoncodeur/cmake-for-complex-projects-part-1-building-a-c-game-engine-from-scratch-for-desktop-and-ied
# Modern CMake magic - precompiled headers for faster builds
# target_precompile_headers(ColumbaEngine PUBLIC src/Engine/stdafx.h)

set( header_files
        include/tdbscan_algo/dummy_logging.h
        include/tdbscan_algo/base_defs.h
        include/tdbscan_algo/absblib.h
        include/tdbscan_algo/common_defs.h
        include/tdbscan_algo/connector.h
        include/tdbscan_algo/connector.hh
        include/tdbscan_algo/tdbcluster.h
        include/tdbscan_algo/tdbcluster.hh
        include/tdbscan_algo/tdbscan_algo.h
        include/tdbscan_algo/tdbscan_algo.hh
)


target_sources(tdbscan_algo
       PUBLIC ${header_files}

)


