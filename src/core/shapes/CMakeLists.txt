file(GLOB Shapes_SRC *.cpp)
add_library(Shapes SHARED ${Shapes_SRC})
install(TARGETS Shapes LIBRARY DESTINATION ${LIBDIR})
set_target_properties(Shapes PROPERTIES SOVERSION ${SOVERSION})
add_dependencies(Shapes EspressoConfig)

