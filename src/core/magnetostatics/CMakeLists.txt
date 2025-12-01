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
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/dipoles.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/dipolar_direct_sum.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/stoner_wohlfarth_thermal.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/dipolar_direct_sum_gpu.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/dlc.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/scafacos_impl.cpp)

if(ESPRESSO_BUILD_WITH_FFTW)
  target_sources(espresso_p3m PRIVATE dp3m.cpp)
  set_source_files_properties(
    dp3m.cpp TARGET_DIRECTORY espresso_p3m
    PROPERTIES INCLUDE_DIRECTORIES ${CMAKE_SOURCE_DIR}/src/core)
endif()
