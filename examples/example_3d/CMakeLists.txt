
add_executable(exec_tdbscan_3d)

target_sources(exec_tdbscan_3d
        PRIVATE
        example_3d.cpp
)

target_link_libraries(exec_tdbscan_3d
        tdbscan_algo
)

target_compile_features(exec_tdbscan_3d PUBLIC cxx_std_17)

set_target_properties(exec_tdbscan_3d PROPERTIES
        CXX_STANDARD_REQUIRED ON
        CXX_EXTENSIONS OFF
)

# Figure out logging: if available use boost logging
if(USE_BOOST_LOGGING)
    target_compile_definitions(
            exec_tdbscan_3d
            PRIVATE
            USE_BOOST_LOGGING=1
    )

    target_link_libraries(exec_tdbscan_3d
            PRIVATE
            Boost::logging
    )
endif()
