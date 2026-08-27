
# --- make available GTest
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

include(FetchContent)
FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG f8d7d77c06936315286eb55f8de22cd23c188571 # release-1.14.0
)
FetchContent_MakeAvailable(googletest)




set(UNIT_TEST ${PROJECT_NAME}_tests)  #local alias

add_executable(${UNIT_TEST} "test_${PROJECT_NAME}.cpp")

target_compile_features(${UNIT_TEST} PRIVATE cxx_std_20)

target_link_libraries(${UNIT_TEST} PRIVATE gtest_main)

add_test(NAME ${UNIT_TEST} COMMAND ${UNIT_TEST})

include_directories("${gtest_SOURCE_DIR}/include")

set(TEST_COMMAND_LINE_OPTIONS ${CMAKE_CTEST_COMMAND} -C $<CONFIGURATION> -R "^${UNIT_TEST}$" --output-on-failures)


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
