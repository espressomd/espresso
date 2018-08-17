set(LBBoundaries_SRC LBBoundary.cpp)
add_library(LBBoundaries SHARED ${LBBoundaries_SRC})
add_dependencies(LBBoundaries EspressoConfig)
target_include_directories(LBBoundaries PRIVATE ${CMAKE_SOURCE_DIR}/src/core ${CMAKE_BINARY_DIR}/src/core)
