file(GLOB Constraints_SRC *.cpp)
add_library(Constraints SHARED ${Constraints_SRC})
add_dependencies(EspressoCore Constraints)

