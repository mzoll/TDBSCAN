#include(FetchContent)
#
#FetchContent_Declare(
#        Catch2
#        GIT_REPOSITORY https://github.com/catchorg/Catch2.git
#        GIT_TAG        v3.4.0 # Replace with the version you need
#)
#
#FetchContent_MakeAvailable(Catch2)

# --- make available GTest
include(FetchContent)
FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest
        GIT_TAG f8d7d77c06936315286eb55f8de22cd23c188571 # release-1.14.0
        #        URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
        #        DOWNLOAD_EXTRACT_TIMESTAMP true
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

# ======================================================
enable_testing()

include(GoogleTest)

add_executable(test_common_clib)
target_link_libraries(test_common_clib
        PRIVATE
        tdbscan_algo
        #Catch2::Catch2
        #Catch2::Catch2WithMain
        GTest::gtest_main
)

target_sources(test_common_clib
        PRIVATE
        ./test_stopwatch.cpp
        ./test_semaphore.cpp
)

gtest_discover_tests(test_common_clib)
#add_test(NAME test_common_clib COMMAND test_common_clib) # Command can be a target

