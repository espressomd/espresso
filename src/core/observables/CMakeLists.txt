file(GLOB Observables_SRC *.?pp)
add_library(Observables SHARED ${Observables_SRC})
add_dependencies(Observables EspressoConfig)
set_target_properties(Observables PROPERTIES MACOSX_RPATH TRUE)
