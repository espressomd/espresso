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

add_library(Espresso_profiler INTERFACE)
add_library(Espresso::profiler ALIAS Espresso_profiler)
target_include_directories(Espresso_profiler INTERFACE "include")

if(WITH_PROFILER)
  find_package(caliper REQUIRED)

  target_link_libraries(Espresso_profiler INTERFACE caliper-mpi)
  target_compile_definitions(Espresso_profiler INTERFACE HAVE_CALIPER)
endif()
