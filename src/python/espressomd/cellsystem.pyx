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
cimport cellsystem
from globals cimport *

cdef class CellSystem(object):
    def set_domain_decomposition(self, use_verlet_lists=True):
        """Activates domain decomposition cell system
        set_domain_decomposition(useVerletList=True)
        """
        if use_verlet_lists:
            dd.use_vList = 1
        else:
            dd.use_vList = 0

        # grid.h::node_grid
        mpi_bcast_cell_structure(CELL_STRUCTURE_DOMDEC)

        # @TODO: gathering should be interface independent
        # return mpi_gather_runtime_errors(interp, TCL_OK)
        return True

    def set_n_square(self, use_verlet_lists=True):
        """Activates the nsquare force calculation
        """
        if use_verlet_lists:
            dd.use_vList = 1
        else:
            dd.use_vList = 0
        mpi_bcast_cell_structure(CELL_STRUCTURE_NSQUARE)
        # @TODO: gathering should be interface independent
        # return mpi_gather_runtime_errors(interp, TCL_OK)
        return True

    def set_layered(self, nLayers=None):
        """set_layered(nLayers=None)
        Set the layerd cell system with nLayers layers"""
        if nLayers:
            if not isinstance(nLayers, int):
                raise ValueError("layer height should be positive")

            if not nLayers > 0:
                raise ValueError("the number of layers has to be >0")

            global n_layers
            n_layers = int(nLayers)
            global determine_n_layers
            determine_n_layers = 0

        if (node_grid[0] != 1 or node_grid[1] != 1):
            node_grid[0] = node_grid[1] = 1
            node_grid[2] = n_nodes
            mpi_err = mpi_bcast_parameter(FIELD_NODEGRID)
        else:
            mpi_err = 0

        if not mpi_err:
            mpi_bcast_cell_structure(CELL_STRUCTURE_LAYERED)

        # @TODO: gathering should be interface independent
        # return mpi_gather_runtime_errors(interp, TCL_OK)

        if mpi_err:
            raise Exception("Broadcasting the node grid failed")
        return True

    def get_state(self):
        s = {}
        if cell_structure.type == CELL_STRUCTURE_LAYERED:
            s["type"] = "layered"
            s["nLayers"] = n_layers
        if cell_structure.type == CELL_STRUCTURE_DOMDEC:
            s["type"] = "domain_decomposition"
            s["use_verlet_lists"] = dd.use_vList
        if cell_structure.type == CELL_STRUCTURE_NSQUARE:
            s["type"] = "nsquare"
            s["use_verlet_lists"] = dd.use_vList

        return s
