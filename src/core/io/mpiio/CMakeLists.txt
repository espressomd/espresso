add_library(mpiio SHARED mpiio.cpp mpiio.hpp)
target_include_directories(mpiio PRIVATE ${CMAKE_SOURCE_DIR}/src/core ${CMAKE_BINARY_DIR}/src/core)
add_dependencies(mpiio EspressoConfig)
