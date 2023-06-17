cmake_minimum_required(VERSION 3.25)
project(vnet_linalg)

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    add_compile_options(-Wall -Wextra -pedantic)
endif()

set(CMAKE_CXX_STANDARD 11)

add_executable(vnet_linalg main.cpp)
add_executable(test_md test_md.cpp)
add_executable(test_kalman test_kalman.cpp)
add_executable(test_static test_static.cpp)