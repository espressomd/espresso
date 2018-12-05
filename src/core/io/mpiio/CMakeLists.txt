add_library(mpiio SHARED mpiio.cpp)
target_include_directories(mpiio PRIVATE ${CMAKE_SOURCE_DIR}/src/core ${CMAKE_BINARY_DIR}/src/core)
target_link_libraries(mpiio PRIVATE EspressoConfig)
