find_package(Boost 1.58 COMPONENTS unit_test_framework REQUIRED)
include(UnitTest)

unit_test(NAME properties_test SRC properties.cpp DEPENDS observables)
unit_test(NAME algorithms_test SRC algorithms.cpp DEPENDS observables)
unit_test(NAME observables_test SRC observables.cpp DEPENDS observables)
