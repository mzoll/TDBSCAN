add_executable(tdbscan_ex)

target_sources(tdbscan_ex
        PRIVATE
        TDBSCAN_3d.cpp)

target_link_libraries(tdbscan_ex
        tdbscan_algo
)


target_compile_features(tdbscan_ex PUBLIC cxx_std_17)

set_target_properties(tdbscan_ex PROPERTIES
        CXX_STANDARD_REQUIRED ON
        CXX_EXTENSIONS OFF
)

#---------------------------------------
add_executable(exec_tdbscan_1d)

target_sources(exec_tdbscan_1d
        PRIVATE
        TDBSCAN_1d.cpp)

target_link_libraries(exec_tdbscan_1d
        tdbscan_algo
        fmt::fmt
)

target_compile_features(exec_tdbscan_1d PUBLIC cxx_std_17)

set_target_properties(exec_tdbscan_1d PROPERTIES
        CXX_STANDARD_REQUIRED ON
        CXX_EXTENSIONS OFF
)



#---------------------------------------
add_executable(exec_tdbscan_3d)

target_sources(exec_tdbscan_3d
        PRIVATE
        TDBSCAN_3d.cpp)

target_link_libraries(exec_tdbscan_3d
        tdbscan_algo
)

target_compile_features(exec_tdbscan_3d PUBLIC cxx_std_17)

set_target_properties(exec_tdbscan_3d PROPERTIES
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


