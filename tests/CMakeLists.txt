# Simple test runner using doctest (single-header test framework)
# https://github.com/doctest/doctest

add_executable(test_solvers test_solvers.cpp)
include_directories(${CMAKE_SOURCE_DIR}/third_party)
target_link_libraries(test_solvers PRIVATE sparse_solver)

enable_testing()

# Register individual test cases with ctest
add_test(NAME Sparse_Matrix_Basic_Operations COMMAND test_solvers --test-case="Sparse Matrix - Basic Operations")
add_test(NAME CG_Solver_Conjugate_Gradient COMMAND test_solvers --test-case="CG Solver - Conjugate Gradient")
add_test(NAME Matrix_Replace_vs_Add COMMAND test_solvers --test-case="Matrix - Replace vs Add")
add_test(NAME Complex_Sparse_Matrix_Basic_Operations COMMAND test_solvers --test-case="Complex Sparse Matrix - Basic Operations")
add_test(NAME COCG_Solver_Conjugate_Orthogonal_CG COMMAND test_solvers --test-case="COCG Solver - Conjugate Orthogonal CG for Symmetric Complex")
