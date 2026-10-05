
add_executable(exec_async)

target_sources(exec_async
        PRIVATE
        example_async.cpp
)

target_link_libraries(exec_async
        tdbscan_algo
)

target_compile_features(exec_async PUBLIC cxx_std_17)

set_target_properties(exec_async PROPERTIES
        CXX_STANDARD_REQUIRED ON
        CXX_EXTENSIONS OFF
)

# Figure out logging: if available use boost logging
if(USE_BOOST_LOGGING)
    target_compile_definitions(
            exec_async
            PRIVATE
            USE_BOOST_LOGGING=1
    )

    target_link_libraries(exec_async
            PRIVATE
            Boost::logging
    )
endif()
