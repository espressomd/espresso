file(GLOB Constraints_SRC *.cpp)
add_library(Constraints SHARED ${Constraints_SRC})
install(TARGETS Constraints LIBRARY DESTINATION ${LIBDIR})
set_target_properties(Constraints PROPERTIES SOVERSION ${SOVERSION})
add_dependencies(Constraints EspressoConfig)

