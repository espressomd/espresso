file(GLOB Shapes_SRC *.cpp)
add_library(Shapes SHARED ${Shapes_SRC})
add_dependencies(Shapes EspressoConfig)
