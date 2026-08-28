add_subdirectory(tdbscan_algo)

add_executable(hello_world
        "../_legacy/old_algorithmic_waste.cxx")

target_sources(hello_world
        PRIVATE
        TDBSCAN.cpp)

target_link_libraries(hello_world
        tdbscan_algo
)

#target_include_libraries()

target_compile_features(hello_world PUBLIC cxx_std_17)

set_target_properties(hello_world PROPERTIES
        CXX_STANDARD_REQUIRED ON
        CXX_EXTENSIONS OFF
)


