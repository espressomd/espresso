include(unit_test)

unit_test(NAME properties_test SRC properties.cpp DEPENDS
          espresso::particle_observables)
unit_test(NAME algorithms_test SRC algorithms.cpp DEPENDS
          espresso::particle_observables)
unit_test(NAME observables_test SRC observables.cpp DEPENDS
          espresso::particle_observables)
