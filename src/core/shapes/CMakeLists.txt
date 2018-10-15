set(Shapes_SRC
    Cylinder.cpp
    Ellipsoid.cpp
    HollowCone.cpp
    Maze.cpp
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
install(TARGETS Shapes LIBRARY DESTINATION ${LIBDIR})
set_target_properties(Shapes PROPERTIES SOVERSION ${SOVERSION})
add_dependencies(Shapes EspressoConfig)
target_include_directories(Shapes PRIVATE ${CMAKE_SOURCE_DIR}/src/core ${CMAKE_BINARY_DIR}/src/core)
