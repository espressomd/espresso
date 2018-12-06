set(LBBoundaries_SRC LBBoundary.cpp)
add_library(LBBoundaries SHARED ${LBBoundaries_SRC})
target_link_libraries(LBBoundaries PRIVATE EspressoConfig Shapes)
target_include_directories(LBBoundaries PRIVATE ${CMAKE_SOURCE_DIR}/src/core ${CMAKE_BINARY_DIR}/src/core)
