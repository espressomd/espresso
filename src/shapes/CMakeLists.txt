add_library(
  espresso_shapes SHARED
  src/HollowConicalFrustum.cpp src/Cylinder.cpp src/Ellipsoid.cpp
  src/Rhomboid.cpp src/Shape.cpp src/SimplePore.cpp src/Slitpore.cpp
  src/Sphere.cpp src/SpheroCylinder.cpp src/Torus.cpp src/Wall.cpp)
add_library(espresso::shapes ALIAS espresso_shapes)

target_link_libraries(espresso_shapes PUBLIC espresso::utils
                      PRIVATE Boost::boost espresso::cpp_flags)
target_include_directories(
  espresso_shapes PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
                         $<INSTALL_INTERFACE:include>)

install(TARGETS espresso_shapes
        LIBRARY DESTINATION ${ESPRESSO_INSTALL_PYTHON}/espressomd)

if(WITH_TESTS)
  add_subdirectory(unit_tests)
endif(WITH_TESTS)
