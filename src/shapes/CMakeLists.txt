set(SOURCE_FILES
    src/HollowConicalFrustum.cpp src/Cylinder.cpp src/Ellipsoid.cpp
    src/Rhomboid.cpp src/SimplePore.cpp src/Slitpore.cpp src/Sphere.cpp
    src/SpheroCylinder.cpp src/Torus.cpp src/Wall.cpp)

add_library(EspressoShapes SHARED ${SOURCE_FILES})
target_link_libraries(EspressoShapes PUBLIC EspressoUtils PRIVATE Boost::boost
                                                                  cxx_interface)
target_include_directories(
  EspressoShapes PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
                        $<INSTALL_INTERFACE:include>)

install(TARGETS EspressoShapes LIBRARY DESTINATION ${PYTHON_INSTDIR}/espressomd)

if(WITH_TESTS)
  add_subdirectory(unit_tests)
endif(WITH_TESTS)
