#
# Copyright (C) 2016-2022 The ESPResSo project
#
# This file is part of ESPResSo.
#
# ESPResSo is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ESPResSo is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

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
