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
from .grid cimport node_grid
from . cimport integrate
from .globals cimport FIELD_SKIN, FIELD_NODEGRID, FIELD_MAXNUMCELLS, FIELD_MINNUMCELLS
from .globals cimport verlet_reuse, skin
from .globals cimport mpi_bcast_parameter
from .cellsystem cimport dd, cell_structure
import numpy as np
from .utils cimport handle_errors
from .utils import is_valid_type

cdef class CellSystem:
    def set_domain_decomposition(self, use_verlet_lists=True,
                                 fully_connected=[False, False, False]):
        """
        Activates domain decomposition cell system.

        Parameters
        ----------
        use_verlet_lists : :obj:`bool`, optional
            Activates or deactivates the usage of Verlet lists
            in the algorithm.

        """

        cell_structure.use_verlet_list = use_verlet_lists
        dd.fully_connected = fully_connected
        # grid.h::node_grid
        mpi_bcast_cell_structure(CELL_STRUCTURE_DOMDEC)

        handle_errors("Error while initializing the cell system.")
        return True

    def set_n_square(self, use_verlet_lists=True):
        """
        Activates the nsquare force calculation.

        Parameters
        ----------
        use_verlet_lists : :obj:`bool`, optional
            Activates or deactivates the usage of the Verlet
            lists for this algorithm.

        """
        cell_structure.use_verlet_list = use_verlet_lists

        mpi_bcast_cell_structure(CELL_STRUCTURE_NSQUARE)
        # @TODO: gathering should be interface independent
        # return mpi_gather_runtime_errors(interp, TCL_OK)
        return True

    def get_state(self):
        s = {"use_verlet_list": cell_structure.use_verlet_list}

        if cell_structure.type == CELL_STRUCTURE_DOMDEC:
            s["type"] = "domain_decomposition"
        if cell_structure.type == CELL_STRUCTURE_NSQUARE:
            s["type"] = "nsquare"

        s["skin"] = skin
        s["verlet_reuse"] = verlet_reuse
        s["n_nodes"] = n_nodes
        s["node_grid"] = np.array([node_grid[0], node_grid[1], node_grid[2]])
        s["cell_grid"] = np.array(
            [dd.cell_grid[0], dd.cell_grid[1], dd.cell_grid[2]])
        s["cell_size"] = np.array(
            [dd.cell_size[0], dd.cell_size[1], dd.cell_size[2]])
        s["fully_connected"] = dd.fully_connected

        return s

    def __getstate__(self):
        s = {"use_verlet_list": cell_structure.use_verlet_list}

        if cell_structure.type == CELL_STRUCTURE_DOMDEC:
            s["type"] = "domain_decomposition"
        if cell_structure.type == CELL_STRUCTURE_NSQUARE:
            s["type"] = "nsquare"

        s["skin"] = skin
        s["node_grid"] = np.array([node_grid[0], node_grid[1], node_grid[2]])
        s["fully_connected"] = dd.fully_connected
        return s

    def __setstate__(self, d):
        use_verlet_lists = None
        for key in d:
            if key == "use_verlet_list":
                use_verlet_lists = d[key]
            elif key == "type":
                if d[key] == "domain_decomposition":
                    self.set_domain_decomposition(
                        use_verlet_lists=use_verlet_lists)
                elif d[key] == "nsquare":
                    self.set_n_square(use_verlet_lists=use_verlet_lists)
        self.skin = d['skin']
        self.node_grid = d['node_grid']

    def get_pairs_(self, distance):
        return mpi_get_pairs(distance)

    def resort(self, global_flag=True):
        """
        Resort the particles in the cellsystem.
        Returns the particle numbers on the nodes
        after the resort.

        Parameters
        ----------
        global_flag : :obj:`bool`
            If true, a global resorting is done, otherwise particles
            are only exchanged between neighboring nodes.

        """

        return mpi_resort_particles(int(global_flag))

    # setter deprecated
    property node_grid:
        """
        Node grid.

        """

        def __set__(self, _node_grid):
            if not np.prod(_node_grid) == n_nodes:
                raise ValueError("Number of available nodes " + str(
                    n_nodes) + " and imposed node grid " + str(_node_grid) + " do not agree.")
            else:
                node_grid[0] = _node_grid[0]
                node_grid[1] = _node_grid[1]
                node_grid[2] = _node_grid[2]
                mpi_err = mpi_bcast_parameter(FIELD_NODEGRID)
                handle_errors("mpi_bcast_parameter for node_grid failed")
                if mpi_err:
                    raise Exception("Broadcasting the node grid failed")

        def __get__(self):
            return np.array([node_grid[0], node_grid[1], node_grid[2]])

    property skin:
        """
        Value of the skin layer expects a floating point number.

        .. note:: Mandatory to set.

        """

        def __set__(self, double _skin):
            if _skin < 0:
                raise ValueError("Skin must be >= 0")
            global skin
            skin = _skin
            mpi_bcast_parameter(FIELD_SKIN)
            integrate.skin_set = True

        def __get__(self):
            return skin

    def tune_skin(self, min_skin=None, max_skin=None, tol=None,
                  int_steps=None, adjust_max_skin=False):
        """
        Tunes the skin by measuring the integration time and bisecting over the
        given range of skins. The best skin is set in the simulation core.

        Parameters
        -----------
        min_skin : :obj:`float`
            Minimum skin to test.
        max_skin : :obj:`float`
            Maximum skin.
        tol : :obj:`float`
            Accuracy in skin to tune to.
        int_steps : :obj:`int`
            Integration steps to time.
        adjust_max_skin : :obj:`bool`, optional
            If ``True``, the value of ``max_skin`` is reduced
            to the maximum permissible skin (in case the passed
            value is too large). Set to ``False`` by default.

        Returns
        -------
        :obj:`float` :
            The :attr:`skin`

        """
        c_tune_skin(min_skin, max_skin, tol, int_steps, adjust_max_skin)
        handle_errors("Error during tune_skin")
        return self.skin
