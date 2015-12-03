
# compile the executable
add_executable(hbt hbt.cpp urqmd_reader.cpp)
target_link_libraries(hbt ${ZLIB_LIBRARIES} ${Boost_LIBRARIES} ${HDF5_CXX_LIBRARIES})
