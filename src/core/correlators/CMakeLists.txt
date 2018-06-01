file(GLOB Correlators_SRC *.?pp)
add_library(Correlators SHARED ${Correlators_SRC})
install(TARGETS Correlators LIBRARY DESTINATION ${LIBDIR})
set_target_properties(Correlators PROPERTIES SOVERSION ${SOVERSION})
add_dependencies(Correlators EspressoConfig)

