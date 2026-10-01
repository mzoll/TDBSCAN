include(FetchContent)

FetchContent_Declare(
        Catch2
        GIT_REPOSITORY github.com/catchorg/Catch2.git
        GIT_TAG        v3.4.0 # Replace with the version you need
)

FetchContent_MakeAvailable(Catch2)
