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

target_sources(
  espresso_core
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/CylindricalLBFluxDensityProfileAtParticlePositions.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/CylindricalLBVelocityProfileAtParticlePositions.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/CylindricalLBVelocityProfile.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/LBVelocityProfile.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/PidObservable.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/RDF.cpp)
