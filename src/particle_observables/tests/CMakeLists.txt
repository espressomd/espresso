include(unit_test)

unit_test(NAME properties_test SRC properties.cpp DEPENDS
          EspressoParticleObservables)
unit_test(NAME algorithms_test SRC algorithms.cpp DEPENDS
          EspressoParticleObservables)
unit_test(NAME observables_test SRC observables.cpp DEPENDS
          EspressoParticleObservables)
