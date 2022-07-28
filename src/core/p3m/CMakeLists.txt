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

if(FFTW3_FOUND)
  target_link_libraries(Espresso_core PUBLIC FFTW3::FFTW3)
endif()

target_sources(
  Espresso_core
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/common.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/TuningAlgorithm.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/send_mesh.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/fft.cpp)
