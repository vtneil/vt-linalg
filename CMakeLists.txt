cmake_minimum_required(VERSION 3.25)
project(vt_linalg)

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    add_compile_options(-Wall -Wextra -Wpedantic)
endif()

set(CMAKE_CXX_STANDARD 11)

include_directories(include)

add_executable(test_md test/test_md.cpp)
add_executable(test_kalman test/test_kalman.cpp)
add_executable(test_static test/test_static.cpp)
add_executable(test_utils test/test_utils.cpp)
add_executable(test_complex test/test_complex.cpp)
add_executable(examples_1d_constant "examples/1-dimension constant rate/main.cpp")
add_executable(examples_3d_tracking "examples/3-dimension tracking/main.cpp")
add_executable(test_compiler test/test_compiler.cpp)
add_executable(main "examples/2-dimension tracking/main.cpp")
add_executable(test_nested_matrix test/test_nested_matrix.cpp)

add_executable(test_tracking test/test_tracking.cpp)
add_executable(test_creation test/test_creation.cpp)