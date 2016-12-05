file(GLOB LBBoundaries_SRC LBBoundary.cpp)
add_library(LBBoundaries SHARED ${LBBoundaries_SRC})
add_dependencies(LBBoundaries EspressoConfig)
set_target_properties(LBBoundaries PROPERTIES MACOSX_RPATH TRUE)
