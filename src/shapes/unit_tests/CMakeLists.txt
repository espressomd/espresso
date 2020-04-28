include(unit_test)
unit_test(NAME Wall_test SRC Wall_test.cpp DEPENDS shapes utils)
unit_test(NAME HollowConicalFrustum_test SRC HollowConicalFrustum_test.cpp
          DEPENDS shapes utils)
unit_test(NAME Union_test SRC Union_test.cpp DEPENDS shapes utils)
unit_test(NAME Ellipsoid_test SRC Ellipsoid_test.cpp DEPENDS shapes utils)
