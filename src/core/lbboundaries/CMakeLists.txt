set(LBBoundaries_SRC LBBoundary.cpp)
add_library(LBBoundaries SHARED ${LBBoundaries_SRC})
add_dependencies(LBBoundaries EspressoConfig)

