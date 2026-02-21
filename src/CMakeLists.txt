add_library(sparse_solver
    Sparse.cpp
    SymSparseProb.cpp
)

target_include_directories(sparse_solver PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_compile_options(sparse_solver PRIVATE
    -Wall -Wextra -O3
)
