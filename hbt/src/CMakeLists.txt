
# build
add_executable(hbt hbt.cpp urqmd_reader.cpp)

# link libraries
target_link_libraries(hbt libhdf5.so libhdf5_cpp.so)
