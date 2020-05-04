set(SOURCE_FILES
    src/HollowConicalFrustum.cpp src/Cylinder.cpp src/Ellipsoid.cpp
    src/Rhomboid.cpp src/SimplePore.cpp src/Slitpore.cpp src/Sphere.cpp
    src/SpheroCylinder.cpp src/Stomatocyte.cpp src/Torus.cpp src/Wall.cpp)

add_library(shapes SHARED ${SOURCE_FILES})
target_link_libraries(shapes PUBLIC utils PRIVATE Boost::boost cxx_interface)
target_include_directories(
  shapes PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
                $<INSTALL_INTERFACE:include>)

install(TARGETS shapes LIBRARY DESTINATION ${PYTHON_INSTDIR}/espressomd)

if(WITH_UNIT_TESTS)
  add_subdirectory(unit_tests)
endif(WITH_UNIT_TESTS)
