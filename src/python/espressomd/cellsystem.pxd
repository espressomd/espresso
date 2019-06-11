#
# Copyright (C) 2013-2018 The ESPResSo project
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

# @TODO: shouldn't these global definitions be used via global_variables?
from __future__ import print_function, absolute_import

from libcpp.vector cimport vector
from libcpp.pair cimport pair

from .utils cimport Vector3i

cdef extern from "communication.hpp":
    void mpi_bcast_cell_structure(int cs)
    int n_nodes
    vector[int] mpi_resort_particles(int global_flag)

cdef extern from "cells.hpp":
    int CELL_STRUCTURE_CURRENT
    int CELL_STRUCTURE_DOMDEC
    int CELL_STRUCTURE_NSQUARE
    int CELL_STRUCTURE_LAYERED

    vector[pair[int, int]] mpi_get_pairs(double distance)

cdef extern from "layered.hpp":
    int determine_n_layers
    int n_layers_ "n_layers"

cdef extern from "tuning.hpp":
    cdef void c_tune_skin "tune_skin" (double min, double max, double tol, int steps)
