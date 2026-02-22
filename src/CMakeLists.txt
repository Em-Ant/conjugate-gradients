add_library(sparse_solver
    sparse_matrix.cpp
    linear_solver_real.cpp
    sparse_matrix_complex.cpp
    linear_solver_complex.cpp
)

target_include_directories(sparse_solver PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_compile_options(sparse_solver PRIVATE
    -Wall -Wextra -O3
)
