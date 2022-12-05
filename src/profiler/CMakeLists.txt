#
# Copyright (C) 2019-2022 The ESPResSo project
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

add_library(espresso_profiler INTERFACE)
add_library(espresso::profiler ALIAS espresso_profiler)
target_include_directories(espresso_profiler INTERFACE "include")

if(ESPRESSO_BUILD_WITH_PROFILER)
  find_package(caliper REQUIRED)

  target_link_libraries(espresso_profiler INTERFACE caliper-mpi)
  target_compile_definitions(espresso_profiler INTERFACE HAVE_CALIPER)
endif()
