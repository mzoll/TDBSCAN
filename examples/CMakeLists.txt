add_executable(tdbscan_ex)

target_sources(tdbscan_ex
        PRIVATE
        TDBSCAN.cpp)

target_link_libraries(tdbscan_ex
        tdbscan_algo
)

#target_include_libraries()

target_compile_features(tdbscan_ex PUBLIC cxx_std_17)

set_target_properties(tdbscan_ex PROPERTIES
        CXX_STANDARD_REQUIRED ON
        CXX_EXTENSIONS OFF
)


