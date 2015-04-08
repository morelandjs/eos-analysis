# cmake version
cmake_minimum_required(VERSION 2.8)

# project name
project(hbt)

# compiler flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++11 -march=native" )

# build
add_subdirectory(src)
