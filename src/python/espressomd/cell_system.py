#
# Copyright (C) 2013-2022 The ESPResSo project
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
from . import utils
from . import particle_data
from .script_interface import ScriptInterfaceHelper, script_interface_register


@script_interface_register
class CellSystem(ScriptInterfaceHelper):
    """
    This class controls the particle decomposition.

    Attributes
    ----------
    decomposition_type : :obj:`str`
        Name of the currently active particle decomposition.
    use_verlet_lists : :obj:`bool`
        Whether to use Verlet lists.
    skin : :obj:`float`
        Verlet list skin.
    node_grid : (3,) array_like of :obj:`int`
        MPI repartition for the regular decomposition cell system.
    max_cut_bonded : :obj:`float`
        Maximal range from bonded interactions.
    max_cut_nonbonded : :obj:`float`
        Maximal range from non-bonded interactions.
    interaction_range : :obj:`float`
        Maximal interaction range from all interactions,
        or -1 when no interactions are active (or their
        cutoff has no impact when only 1 MPI rank is used).

    Methods
    -------
    resort()
        Resort the particles in the cell system.

        Parameters
        ----------
        global_flag : :obj:`bool`, optional
            If true (default), a global resorting is done, otherwise particles
            are only exchanged between neighboring nodes.

        Returns
        -------
        (N,) array_like of :obj:`int`
            The number of particles per node.

    tune_skin()
        Tune the skin by measuring the integration time and bisecting over the
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
            value is too large). Defaults to ``False``.

        Returns
        -------
        :obj:`float` :
            The :attr:`skin`

    get_state()
        Get the current state of the cell system.

        Returns
        -------
        :obj:`dict` :
            The cell system state.

    """
    _so_name = "CellSystem::CellSystem"
    _so_creation_policy = "GLOBAL"
    _so_bind_methods = ("get_state", "tune_skin", "resort")

    def set_regular_decomposition(self, **kwargs):
        """
        Activate the regular decomposition cell system.

        Parameters
        ----------
        use_verlet_lists : :obj:`bool`, optional
            Activates or deactivates the usage of Verlet lists.
            Defaults to ``True``.
        fully_connected_boundary : :obj:`dict`, optional
            If set, connects all cells on a given boundary along the given direction.
            Example: ``{"boundary": "z", "direction": "x"}`` connects all
            cells on the boundary normal to the z-direction along the x-axis.
            This corresponds to z-axis as shear plane normal and x-axis as
            shear direction in Lees-Edwards boundary conditions.

        """
        self.call_method("initialize", name="regular_decomposition", **kwargs)

    def set_n_square(self, **kwargs):
        """
        Activate the N-square cell system.

        Parameters
        ----------
        use_verlet_lists : :obj:`bool`, optional
            Activates or deactivates the usage of Verlet lists.
            Defaults to ``True``.

        """
        self.call_method("initialize", name="n_square", **kwargs)

    def set_hybrid_decomposition(self, **kwargs):
        """
        Activate the hybrid domain decomposition.

        Parameters
        ----------
        cutoff_regular : :obj:`float`
            Maximum cutoff to consider in regular decomposition.
            Should be as low as the system permits.
        n_square_types : list of :obj:`int`, optional
            Particles types that should be handled in the N-square cell system.
            Defaults to an empty list.
        use_verlet_lists : :obj:`bool`, optional
            Activates or deactivates the usage of Verlet lists.
            Defaults to ``True``.

        """
        self.call_method("initialize", name="hybrid_decomposition", **kwargs)

    def get_pairs(self, distance, types='all'):
        """
        Get pairs of particles closer than threshold value.

        Parameters
        ----------
        distance : :obj:`float`
            Pairs of particles closer than ``distance`` are found.
        types : list of :obj:`int` or ``'all'``, optional
            Restrict the pair search to the specified types.
            Defaults to ``'all'``, in which case all particles are considered.

        Returns
        -------
        list of tuples of :obj:`int`
            The particle pairs identified by their index

        Raises
        ------
        Exception
            If the pair search distance is greater than the cell size
        """
        pairs = self.call_method("get_pairs", distance=distance, types=types)
        return [tuple(pair) for pair in pairs]

    def get_neighbors(self, particle, distance):
        """
        Get neighbors of a given particle up to a certain distance.

        The choice of :ref:`cell systems <Cell systems>` has an impact on
        how far the algorithm can detect particle pairs:

        * N-square: no restriction on the search distance, no double counting
          if search distance is larger than the box size
        * regular decomposition: the search distance is bounded by half
          the local cell geometry
        * hybrid decomposition: not supported

        Parameters
        ----------
        particle : :class:`~espressomd.particle_data.ParticleHandle`
        distance : :obj:`float`
            Pairs of particles closer than ``distance`` are found.

        Returns
        -------
        (N,) array_like of :obj:`int`
            The list of neighbor particles surrounding the particle

        """
        utils.check_type_or_throw_except(
            particle, 1, particle_data.ParticleHandle, "Parameter 'particle' must be a particle")
        utils.check_type_or_throw_except(
            distance, 1, float, "Parameter 'distance' must be a float")
        return self.call_method(
            "get_neighbors", distance=distance, pid=particle.id)

    def non_bonded_loop_trace(self):
        pairs = self.call_method("non_bonded_loop_trace")
        res = []
        for id1, id2, pos1, pos2, vec21, node in pairs:
            pos1 = utils.array_locked(pos1)
            pos2 = utils.array_locked(pos2)
            vec21 = utils.array_locked(vec21)
            res.append((id1, id2, pos1, pos2, vec21, node))
        return res
