#
# Copyright (C) 2010-2022 The ESPResSo project
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

from .script_interface import ScriptInterfaceHelper, script_interface_register


@script_interface_register
class CollisionDetection(ScriptInterfaceHelper):

    """
    Interface to the collision detection / dynamic binding.

    See :ref:`Creating bonds when particles collide` for detailed instructions.

    This class should not be instantiated by the user. Instead, use
    the :attr:`~espressomd.system.System.collision_detection` attribute
    of the system class to access the collision detection.

    Use method :meth:`~espressomd.collision_detection.CollisionDetection.set_params`
    to change the parameters of the collision detection.

    """

    _so_name = "CollisionDetection::CollisionDetection"
    _so_features = ("COLLISION_DETECTION",)

    # Do not allow setting of individual attributes
    def __setattr__(self, *args, **kwargs):
        raise Exception(
            "Please set all parameters at once via collision_detection.set_params()")

    def set_params(self, **kwargs):
        """
        Set the parameters for the collision detection

        See :ref:`Creating bonds when particles collide` for detailed instructions.


        Parameters
        ----------
        mode : :obj:`str`, {"off", "bind_centers", "bind_at_point_of_collision", "glue_to_surface"}
            Collision detection mode

        distance : :obj:`float`
            Distance below which a pair of particles is considered in the
            collision detection

        bond_centers : :obj:`espressomd.interactions.BondedInteraction`
            Bond to add between the colliding particles

        bond_vs : :obj:`espressomd.interactions.BondedInteraction`
            Bond to add between virtual sites (for modes using virtual sites)

        part_type_vs : :obj:`int`
            Particle type of the virtual sites being created on collision
            (virtual sites based modes)

        part_type_to_be_glued : :obj:`int`
            particle type for ``"glue_to_surface"`` mode. See user guide.

        part_type_to_attach_vs_to : :obj:`int`
            particle type for ``"glue_to_surface"`` mode. See user guide.

        part_type_after_glueing : :obj:`int`
            particle type for ``"glue_to_surface"`` mode. See user guide.

        distance_glued_particle_to_vs : :obj:`float`
            Distance for ``"glue_to_surface"`` mode. See user guide.

        """

        if "mode" not in kwargs:
            raise ValueError(
                "Collision mode must be specified via the 'mode' argument")
        self.call_method("set_params", **kwargs)

    def get_parameter(self, name):
        value = super().get_parameter(name)
        if name in ["bond_centers", "bond_vs"]:
            if value == -1:  # Not defined
                value = None
            else:
                value = self.call_method("get_bond_by_id", bond_id=value)
        return value

    def get_params(self):
        """Returns the parameters of the collision detection as dict.

        """
        params = {}
        mode = super().get_parameter("mode")
        for name in self.call_method("params_for_mode", mode=mode):
            params[name] = self.get_parameter(name)
        return params
