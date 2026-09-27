
add_executable(exec_tdbscan_1d)

target_sources(exec_tdbscan_1d
        PRIVATE
        example_1d.cpp
)

target_include_directories(exec_tdbscan_1d
        PRIVATE
        ./include
)

target_link_libraries(exec_tdbscan_1d
        tdbscan_algo
        fmt::fmt
)

target_compile_features(exec_tdbscan_1d PUBLIC cxx_std_17)

set_target_properties(exec_tdbscan_1d PROPERTIES
        CXX_STANDARD_REQUIRED ON
        CXX_EXTENSIONS OFF
)
