set(LBBoundaries_SRC LBBoundary.cpp)
add_library(LBBoundaries SHARED ${LBBoundaries_SRC})
target_link_libraries(LBBoundaries PUBLIC EspressoConfig EspressoCore Shapes)
