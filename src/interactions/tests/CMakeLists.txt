include(unit_test)

unit_test(NAME lennard_jones_test SRC lennard_jones.test.cpp DEPENDS
          EspressoInteractions)
unit_test(NAME linear_test SRC linear.test.cpp DEPENDS EspressoInteractions)
unit_test(NAME central_potential_test SRC central_potential.test.cpp DEPENDS
          EspressoInteractions)
unit_test(NAME container_test SRC registry.test.cpp DEPENDS
          EspressoInteractions)
