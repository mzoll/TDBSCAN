add_library(tdbscan_algo
        src/include_all.cpp
)

target_include_directories(tdbscan_algo
        PUBLIC include
        PRIVATE src
)

# /// IDEA ///
# Modern CMake magic - precompiled headers for faster builds
# target_precompile_headers(<target> PUBLIC <header-files>)

list(APPEND header_files
        include/tdbscan_algo/auxilary/auxil.h
        include/tdbscan_algo/auxilary/trivial_logging.h
        include/tdbscan_algo/base_defs.h
        include/tdbscan_algo/absblib.h
        include/tdbscan_algo/common_defs.h
        #include/tdbscan_algo/common_defs.hh
        include/tdbscan_algo/common_blibs.h
        #include/tdbscan_algo/common_defs.hh
        include/tdbscan_algo/connector.h
        #include/tdbscan_algo/connector.hh
        include/tdbscan_algo/tdbcluster.h
        #include/tdbscan_algo/tdbcluster.hh
        include/tdbscan_algo/tdbscan_algo.h
        #include/tdbscan_algo/tdbscan_algo.hh
)

target_sources(tdbscan_algo
       PUBLIC ${header_files}
       #src/include_all.cpp
)

# Figure out logging: if available use boost logging
if(USE_BOOST_LOGGING)
    target_compile_definitions(
        tdbscan_algo
        INTERFACE
            USE_BOOST_LOGGING=1
    )

    target_link_libraries(tdbscan_algo
        PUBLIC
        Boost::logging
    )
endif()
