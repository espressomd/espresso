# Copyright (C) 2010-2020 The ESPResSo project
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

from libcpp.string cimport string
from libcpp.vector cimport vector

cdef extern from "communication.hpp":
    cdef void mpi_bcast_cell_structure(int cs)
    cdef void mpi_bcast_generic_dd_grid(string)
    int CELL_STRUCTURE_GENERIC_DD

cdef extern from "generic-dd/generic_dd.hpp" namespace "generic_dd":
    cdef vector[string] librepa_supported_grid_types()
