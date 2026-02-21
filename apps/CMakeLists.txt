# Simple test runner that mimics original main.cpp functionality
add_executable(test_runner test_runner.cpp)
target_link_libraries(test_runner PRIVATE sparse_solver)
