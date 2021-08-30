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
import numpy as np
from libcpp.cast cimport dynamic_cast
from .grid cimport node_grid
from .grid cimport box_geo
from . cimport integrate
from libcpp.vector cimport vector
from .cellsystem cimport cell_structure
from .cellsystem cimport get_verlet_reuse
from .cellsystem cimport mpi_set_skin, skin
from .utils import handle_errors
from .utils cimport Vector3i
from .utils cimport check_type_or_throw_except, make_array_locked

cdef class CellSystem:
    def set_domain_decomposition(self, use_verlet_lists=True):
        """
        Activates domain decomposition cell system.

        Parameters
        ----------
        use_verlet_lists : :obj:`bool`, optional
            Activates or deactivates the usage of Verlet lists
            in the algorithm.

        """
        mpi_set_use_verlet_lists(use_verlet_lists)
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
        mpi_set_use_verlet_lists(use_verlet_lists)
        mpi_bcast_cell_structure(CELL_STRUCTURE_NSQUARE)

        return True

    def get_state(self):
        s = self.__getstate__()

        if cell_structure.decomposition_type() == CELL_STRUCTURE_DOMDEC:
            dd = get_domain_decomposition()
            s["cell_grid"] = np.array(
                [dd.cell_grid[0], dd.cell_grid[1], dd.cell_grid[2]])
            s["cell_size"] = np.array(
                [dd.cell_size[0], dd.cell_size[1], dd.cell_size[2]])

        s["verlet_reuse"] = get_verlet_reuse()
        s["n_nodes"] = n_nodes

        return s

    def __getstate__(self):
        s = {"use_verlet_list": cell_structure.use_verlet_list}

        if cell_structure.decomposition_type() == CELL_STRUCTURE_DOMDEC:
            s["type"] = "domain_decomposition"
        if cell_structure.decomposition_type() == CELL_STRUCTURE_NSQUARE:
            s["type"] = "nsquare"

        s["skin"] = skin
        s["node_grid"] = np.array([node_grid[0], node_grid[1], node_grid[2]])
        return s

    def __setstate__(self, d):
        self.skin = d['skin']
        self.node_grid = d['node_grid']
        if 'type' in d:
            if d['type'] == "domain_decomposition":
                self.set_domain_decomposition(
                    use_verlet_lists=d['use_verlet_list'])
            elif d['type'] == "nsquare":
                self.set_n_square(use_verlet_lists=d['use_verlet_list'])

    def get_pairs(self, distance, types='all'):
        """
        Get pairs of particles closer than threshold value

        Parameters
        ----------
        distance : :obj:`float`
            Pairs of particles closer than ``distance`` are found.
        types (optional) : list of :obj:`int` or ``'all'``
            Restrict the pair search to the specified types. Defaults to ``'all'``, in which case all particles are considered.

        Returns
        -------
        list of tuples of :obj:`int`
            The particle pairs identified by their index

        Raises
        ------
        Exception
            If the pair search distance is greater than the cell size
        """
        pairs = None
        if types == 'all':
            pairs = mpi_get_pairs(distance)
        else:
            pairs = self._get_pairs_of_types(distance, types)
        handle_errors("")
        return pairs

    def _get_pairs_of_types(self, distance, types):
        """
        This function needs to be separated from ``self.get_pairs()`` because ``cdef``
        cannot be inside an ``else`` block
        """
        if hasattr(types, "__getitem__"):
            if len(types) == 1:
                check_type_or_throw_except(
                    types[0], 1, int, 'types must be a list of int')
            else:
                check_type_or_throw_except(
                    types, len(types), int, 'types must be a list of int')
        else:
            raise ValueError('types must be iterable')

        cdef vector[int] types_c
        for type in types:
            types_c.push_back(type)
        return mpi_get_pairs_of_types(distance, types_c)

    def non_bonded_loop_trace(self):
        cdef vector[PairInfo] pairs = mpi_non_bonded_loop_trace()
        cdef PairInfo pair
        res = []
        for pair in pairs:
            res.append(
                (pair.id1, pair.id2,
                 make_array_locked(pair.pos1),
                 make_array_locked(pair.pos2),
                 make_array_locked(pair.vec21), pair.node))
        return res

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

    property node_grid:
        """
        Node grid.

        """

        def __set__(self, _node_grid):
            cdef Vector3i node_grid
            if not np.prod(_node_grid) == n_nodes:
                raise ValueError(
                    f"Number of available nodes {n_nodes} and imposed node grid {_node_grid} do not agree.")
            else:
                for i in range(3):
                    node_grid[i] = _node_grid[i]
                mpi_set_node_grid(node_grid)

        def __get__(self):
            return np.array([node_grid[0], node_grid[1], node_grid[2]])

    property max_cut_nonbonded:
        def __get__(self):
            return maximal_cutoff_nonbonded()

    property max_cut_bonded:
        def __get__(self):
            return maximal_cutoff_bonded()

    property skin:
        """
        Value of the skin layer expects a floating point number.

        .. note:: Mandatory to set.

        """

        def __set__(self, double _skin):
            if _skin < 0:
                raise ValueError("Skin must be >= 0")
            mpi_set_skin(_skin)

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
