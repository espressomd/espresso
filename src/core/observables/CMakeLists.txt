file(GLOB Observables_SRC *.?pp)
add_library(Observables SHARED ${Observables_SRC})
install(TARGETS Observables LIBRARY DESTINATION ${LIBDIR})
set_target_properties(Observables PROPERTIES SOVERSION ${SOVERSION})
add_dependencies(Observables EspressoConfig)

