#
# Copyright (C) 2018-2022 The ESPResSo project
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
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/bmhtf-nacl.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/buckingham.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/gaussian.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/gay_berne.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/hat.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/hertzian.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/ljcos2.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/ljcos.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/lj.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/ljgen.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/morse.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/nonbonded_interaction_data.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/soft_sphere.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/smooth_step.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/wca.cpp)
