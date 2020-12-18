#
# Copyright (C) 2013-2019 The ESPResSo project
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

from libcpp cimport bool
from libcpp.vector cimport vector
from libcpp.pair cimport pair
from .utils cimport Vector3i

cdef extern from "communication.hpp":
    int n_nodes

cdef extern from "cells.hpp":
    int CELL_STRUCTURE_DOMDEC
    int CELL_STRUCTURE_NSQUARE

    ctypedef struct CellStructure:
        int decomposition_type()
        bool use_verlet_list

    CellStructure cell_structure

    const DomainDecomposition * get_domain_decomposition()

    vector[pair[int, int]] mpi_get_pairs(double distance) except +
    vector[pair[int, int]] mpi_get_pairs_of_types(double distance, vector[int] types) except +
    vector[int] mpi_resort_particles(int global_flag)
    void mpi_bcast_cell_structure(int cs)
    void mpi_set_use_verlet_lists(bool use_verlet_lists)

cdef extern from "tuning.hpp":
    cdef void c_tune_skin "tune_skin" (double min_skin, double max_skin, double tol, int int_steps, bool adjust_max_skin)

cdef extern from "DomainDecomposition.hpp":
    cppclass  DomainDecomposition:
        Vector3i cell_grid
        double cell_size[3]
