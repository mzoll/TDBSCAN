add_test(
        NAME MyAppWithTestFlag
        COMMAND MyApp --test
)


target_sources(Testsuite
    PRIVATE
    test/pure_algo/test_absblib.cpp
)
