#
# Copyright (C) 2013,2014,2015,2016 The ESPResSo project
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
cdef extern from "tcl.h":
    # @TODO: this struct probably should be defined at some common module
    cdef struct Tcl_Interp:
        char * result
        int errorLine

# @TODO: shouldn't these global definitions be used via global_variables?
cdef extern from "global.hpp":
    int FIELD_NODEGRID

cdef extern from "communication.hpp":
    int mpi_bcast_parameter(int p)
    int mpi_gather_runtime_errors(Tcl_Interp * interp, int ret_state)
    void mpi_bcast_cell_structure(int cs)
    int n_nodes

cdef extern from "cells.hpp":
    int CELL_STRUCTURE_CURRENT
    int CELL_STRUCTURE_DOMDEC
    int CELL_STRUCTURE_NSQUARE
    int CELL_STRUCTURE_LAYERED

cdef extern from "layered.hpp":
    int determine_n_layers
    int n_layers
    int determine_n_layers

cdef extern from "grid.hpp":
    int node_grid[3]
