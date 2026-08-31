
# --- make available GTest
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

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

add_executable(test_common_defs
        tdbscan_algo/test_common_defs.cpp
)
target_link_libraries(test_common_defs
        PRIVATE tdbscan_algo
        GTest::gtest_main
)
#-----------------------------
add_executable(test_aux
        tdbscan_algo/test_aux.cpp
)
target_link_libraries(test_aux
        PRIVATE tdbscan_algo
        GTest::gtest_main
)
#-----------------------------
add_executable(test_connector
        tdbscan_algo/test_connector.cpp
)
target_link_libraries(test_connector
        PRIVATE tdbscan_algo
        GTest::gtest_main
)
#-----------------------------
add_executable(test_absblib
        tdbscan_algo/test_absblib.cpp
)
target_link_libraries(test_absblib
        PRIVATE tdbscan_algo
        GTest::gtest_main
)


gtest_discover_tests(test_common_defs)
gtest_discover_tests(test_connector)
gtest_discover_tests(test_absblib)

# --- adding test through cmake magic
## define tests
#add_test(
#        NAME test_basedefs
#        COMMAND $<TARGET_FILE:test_basedefs>
#)




# ======================================================
#
#
#set(UNIT_TEST ${PROJECT_NAME}_tests)  #local alias
#
#add_executable(${UNIT_TEST} "test_${PROJECT_NAME}.cpp")
#
#target_compile_features(${UNIT_TEST} PRIVATE cxx_std_20)
#
#target_link_libraries(${UNIT_TEST} PRIVATE gtest_main)
#
#add_test(NAME ${UNIT_TEST} COMMAND ${UNIT_TEST})
#
#include_directories("${gtest_SOURCE_DIR}/include")
#
#set(TEST_COMMAND_LINE_OPTIONS ${CMAKE_CTEST_COMMAND} -C $<CONFIGURATION> -R "^${UNIT_TEST}$" --output-on-failures)
#

# ======================================================

#
#add_test(
#        NAME MyAppWithTestFlag
#        COMMAND MyApp --test
#)
#
#
#target_sources(Testsuite
#    PRIVATE
#    test/tdbscan_algo/test_absblib.cpp
#)
