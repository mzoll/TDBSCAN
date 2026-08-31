add_library(tdbscan_algo
        src/include_all.cpp
        )

target_include_directories(tdbscan_algo
        PUBLIC include
        PRIVATE src
)

target_link_libraries(tdbscan_algo

)

# /// IDEA ///
# https://dev.to/pigeoncodeur/cmake-for-complex-projects-part-1-building-a-c-game-engine-from-scratch-for-desktop-and-ied
# Modern CMake magic - precompiled headers for faster builds
# target_precompile_headers(ColumbaEngine PUBLIC src/Engine/stdafx.h)

set( header_files
        tdbscan_algo/dummy_logging.h
        tdbscan_algo/aux.h
        tdbscan_algo/base_defs.h
        tdbscan_algo/absblib.h
        tdbscan_algo/common_defs.h
        tdbscan_algo/tdbcluster.h
        tdbscan_algo/tdbcluster.hh
        tdbscan_algo/connector.h
        tdbscan_algo/connector.hh
        tdbscan_algo/tdbscan_algo.h
        tdbscan_algo/tdbscan_algo.hh
)


target_sources(tdbscan_algo
       PUBLIC ${header_files}

)


