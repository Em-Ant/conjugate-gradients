cmake_minimum_required(VERSION 3.10)
project(ConjugateGradients LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()

if(CMAKE_VERSION VERSION_LESS "3.24")
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
else()
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
endif()

enable_testing()

add_subdirectory(src)
add_subdirectory(tests)
add_subdirectory(apps)
