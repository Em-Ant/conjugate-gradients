add_library(sparse_solver INTERFACE)

target_include_directories(sparse_solver INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_compile_options(sparse_solver INTERFACE
    -Wall -Wextra -O3
)
