include(unit_test)

unit_test(NAME properties_test SRC properties.cpp DEPENDS genobs)
unit_test(NAME algorithms_test SRC algorithms.cpp DEPENDS genobs)
unit_test(NAME observables_test SRC observables.cpp DEPENDS genobs)
