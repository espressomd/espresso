set(Shapes_SRC
        HollowConicalFrustum.cpp
    Cylinder.cpp
    Ellipsoid.cpp
    Rhomboid.cpp
    SimplePore.cpp
    Slitpore.cpp
    Sphere.cpp
    SpheroCylinder.cpp
    Stomatocyte.cpp
    Torus.cpp
    Wall.cpp
)
add_library(Shapes SHARED ${Shapes_SRC})
install(TARGETS Shapes LIBRARY DESTINATION ${PYTHON_INSTDIR}/espressomd)
target_link_libraries(Shapes PUBLIC EspressoConfig EspressoCore)

if(WITH_UNIT_TESTS)
  add_subdirectory(unit_tests)
endif(WITH_UNIT_TESTS)

