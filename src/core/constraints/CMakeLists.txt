file(GLOB Constraints_SRC Constraint.cpp)
add_library(Constraints SHARED ${Constraints_SRC})
add_dependencies(Constraints EspressoConfig)
set_target_properties(Constraints PROPERTIES MACOSX_RPATH TRUE)
