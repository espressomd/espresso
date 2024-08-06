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
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/coulomb.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/elc.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/icc.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/mmm1d.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/mmm-modpsi.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/p3m.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/scafacos_impl.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/specfunc.cpp)
