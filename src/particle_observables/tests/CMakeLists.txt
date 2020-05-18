include(unit_test)

unit_test(NAME properties_test SRC properties.cpp DEPENDS particle_observables)
unit_test(NAME algorithms_test SRC algorithms.cpp DEPENDS particle_observables)
unit_test(NAME observables_test SRC observables.cpp DEPENDS
          particle_observables)
