#---------------------------------------
add_subdirectory(example_1d)

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


