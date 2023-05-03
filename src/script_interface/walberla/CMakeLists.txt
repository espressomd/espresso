#
# Copyright (C) 2021-2023 The ESPResSo project
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

if(ESPRESSO_BUILD_WITH_WALBERLA)
  target_sources(
    espresso_script_interface
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/initialize.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/LBFluidNode.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/LBFluidSlice.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/LBFluid.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/EKSpecies.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/EKSpeciesNode.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/EKSpeciesSlice.cpp)
  target_link_libraries(espresso_script_interface PRIVATE espresso::walberla)
endif(ESPRESSO_BUILD_WITH_WALBERLA)
