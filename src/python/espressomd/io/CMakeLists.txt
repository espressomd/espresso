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

configure_file(mpiio.py mpiio.py COPYONLY)
configure_file(vtk.py vtk.py COPYONLY)
add_subdirectory(writer)
set(cython_AUX ${cython_AUX}
               "${CMAKE_SOURCE_DIR}/src/python/espressomd/io/__init__.py"
    CACHE INTERNAL "cython_AUX" FORCE)
