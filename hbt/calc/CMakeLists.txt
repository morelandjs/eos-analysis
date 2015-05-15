# cmake version
cmake_minimum_required(VERSION 2.8)

# project name
project(hbt)

# find zlib
find_package(ZLIB REQUIRED)

# find boost
find_package(Boost REQUIRED COMPONENTS filesystem system iostreams)

# find hdf5
find_package(HDF5 REQUIRED COMPONENTS CXX)

# compiler flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -static -std=c++11 -march=corei7" )

# add the src directory
add_subdirectory(src)
