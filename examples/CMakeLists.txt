#---------------------------------------
add_subdirectory(example_1d)

#---------------------------------------
add_subdirectory(example_3d)

add_subdirectory(example_async)

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


