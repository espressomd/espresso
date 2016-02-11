
file(GLOB ImmersedBoundary_SRC *.cpp)
add_library(ImmersedBoundary SHARED ${ImmersedBoundary_SRC})
set_target_properties(ImmersedBoundary PROPERTIES MACOSX_RPATH TRUE)

set_target_properties(ImmersedBoundary PROPERTIES LINKER_LANGUAGE CXX)
