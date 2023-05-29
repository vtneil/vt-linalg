cmake_minimum_required(VERSION 3.25)
project(vnet_linalg)

set(CMAKE_CXX_STANDARD 11)

add_executable(vnet_linalg main.cpp)

add_executable(mm_ijk mm_ijk.cpp)