add_library(mpiio SHARED mpiio.cpp)
target_include_directories(mpiio PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(mpiio PRIVATE EspressoConfig EspressoCore MPI::MPI_CXX)
