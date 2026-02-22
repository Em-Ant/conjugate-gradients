cmake_minimum_required(VERSION 3.10)
project(ConjugateGradients LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Default to Debug build for IntelliSense compatibility
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "Build type" FORCE)
endif()

# Generate compile_commands.json for IDE integration
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Include third_party for all targets (doctest test framework)
include_directories(${CMAKE_SOURCE_DIR}/third_party)

enable_testing()
include(CTest)

add_subdirectory(src)
add_subdirectory(tests)
