include(FetchContent)

FetchContent_Declare(
    catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG v3.4.0
)

FetchContent_MakeAvailable(catch2)

include(CTest)
include(Catch)

add_executable(test_solvers test_solvers.cpp)
target_link_libraries(test_solvers PRIVATE sparse_solver Catch2::Catch2WithMain)

catch_discover_tests(test_solvers)
