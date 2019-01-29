set(Shapes_SRC
    Cylinder.cpp
    Ellipsoid.cpp
    HollowCone.cpp
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
set_target_properties(Shapes PROPERTIES SOVERSION ${SOVERSION})
target_link_libraries(Shapes PUBLIC EspressoConfig EspressoCore)
