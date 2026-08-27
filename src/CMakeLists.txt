include(cmake/StandardProjectSettings.cmake)

add_subdirectory(pure_algo)

add_executable(hello_world)

target_sources(hello_world
        PRIVATE
        src/TDBSCAN.cpp)

target_compile_features(hello_world PUBLIC cxx_std_17)

set_target_properties(hello_world PROPERTIES
        CXX_STANDARD_REQUIRED ON
        CXX_EXTENSIONS OFF
)


