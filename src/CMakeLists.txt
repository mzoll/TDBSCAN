add_subdirectory(tdbscan_algo)

add_executable(tdbscan)

target_sources(tdbscan
        PRIVATE
        TDBSCAN.cpp)

target_link_libraries(tdbscan
        tdbscan_algo
)

#target_include_libraries()

target_compile_features(tdbscan PUBLIC cxx_std_17)

set_target_properties(tdbscan PROPERTIES
        CXX_STANDARD_REQUIRED ON
        CXX_EXTENSIONS OFF
)


