add_executable(tdbscan_ex)

target_sources(tdbscan_ex
        PRIVATE
        TDBSCAN.cpp)

target_link_libraries(tdbscan_ex
        tdbscan_algo
)


target_compile_features(tdbscan_ex PUBLIC cxx_std_17)

set_target_properties(tdbscan_ex PROPERTIES
        CXX_STANDARD_REQUIRED ON
        CXX_EXTENSIONS OFF
)

#---------------------------------------
add_executable(tdbscan_ex1d)

target_sources(tdbscan_ex1d
        PRIVATE
        TDBSCAN_1d.cpp)

target_link_libraries(tdbscan_ex1d
        tdbscan_algo
)

target_compile_features(tdbscan_ex1d PUBLIC cxx_std_17)

set_target_properties(tdbscan_ex1d PROPERTIES
        CXX_STANDARD_REQUIRED ON
        CXX_EXTENSIONS OFF
)


#---------------------------------------
add_executable(sandbox)

target_sources(sandbox
        PRIVATE
        sandbox.cpp)

target_link_libraries(sandbox
        tdbscan_algo
)

target_compile_features(sandbox PUBLIC cxx_std_17)

set_target_properties(sandbox PROPERTIES
        CXX_STANDARD_REQUIRED ON
        CXX_EXTENSIONS OFF
)


