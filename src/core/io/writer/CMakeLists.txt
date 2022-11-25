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

if(ESPRESSO_BUILD_WITH_HDF5)
  target_link_libraries(espresso_core PUBLIC ${HDF5_LIBRARIES}
                                             Boost::filesystem h5xx)
  target_include_directories(espresso_core PUBLIC ${CMAKE_CURRRENT_SOURCE_DIR}
                                                  ${HDF5_INCLUDE_DIRS})

  target_compile_definitions(espresso_core PUBLIC H5XX_USE_MPI)

  target_sources(
    espresso_core PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/h5md_core.cpp
                          ${CMAKE_CURRENT_SOURCE_DIR}/h5md_specification.cpp)
endif()
