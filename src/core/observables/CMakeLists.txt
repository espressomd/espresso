file(GLOB Observables_SRC *.?pp)
add_library(Observables SHARED ${Observables_SRC})
add_dependencies(Observables EspressoConfig)

