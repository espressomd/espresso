find_package(Boost 1.58 COMPONENTS unit_test_framework REQUIRED)
include(UnitTest)

unit_test(NAME traits_test SRC traits.cpp DEPENDS traits)
unit_test(NAME algorithms_test SRC algorithms.cpp DEPENDS algorithms)
unit_test(NAME observables_test SRC observables.cpp DEPENDS observables traits)
