add_library(
  Espresso_shapes SHARED
  src/HollowConicalFrustum.cpp src/Cylinder.cpp src/Ellipsoid.cpp
  src/Rhomboid.cpp src/SimplePore.cpp src/Slitpore.cpp src/Sphere.cpp
  src/SpheroCylinder.cpp src/Torus.cpp src/Wall.cpp)
add_library(Espresso::shapes ALIAS Espresso_shapes)

target_link_libraries(Espresso_shapes PUBLIC Espresso::utils
                      PRIVATE Boost::boost Espresso::cpp_flags)
target_include_directories(
  Espresso_shapes PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
                         $<INSTALL_INTERFACE:include>)

install(TARGETS Espresso_shapes
        LIBRARY DESTINATION ${PYTHON_INSTDIR}/espressomd)

if(WITH_TESTS)
  add_subdirectory(unit_tests)
endif(WITH_TESTS)
