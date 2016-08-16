file(GLOB Constraints_SRC Constraint.cpp)
add_library(Constraints SHARED ${Constraints_SRC})
add_dependencies(Constraints EspressoConfig)

