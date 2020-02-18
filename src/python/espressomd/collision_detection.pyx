# Copyright (C) 2010-2019 The ESPResSo project
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
from .script_interface import ScriptInterfaceHelper, script_interface_register
from .utils import to_str
from .utils cimport handle_errors
from .interactions import BondedInteraction, BondedInteractions


cdef extern from "collision.hpp":
    const int COLLISION_MODE_OFF
    const int COLLISION_MODE_BOND
    const int COLLISION_MODE_VS
    const int COLLISION_MODE_GLUE_TO_SURF
    const int COLLISION_MODE_BIND_THREE_PARTICLES


@script_interface_register
class CollisionDetection(ScriptInterfaceHelper):

    """
    Interface to the collision detection / dynamic binding.

    See :ref:`Creating bonds when particles collide` for detailed instructions.

    This class should not be instantiated by the user. Instead, use
    the :attr:`espressomd.system.System.collision_detection` attribute
    of the system class to access the collision detection.

    Use :meth:`espressomd.collision_detection.CollisionDetection.set_params`
    to change the parameters of the collision detection.

    """

    _so_name = "CollisionDetection::CollisionDetection"

    def __init__(self, *args, **kwargs):
        super().__init__()
        self.set_params(**kwargs)

    def validate(self):
        """Validates the parameters of the collision detection.

        This is called automatically on parameter change

        """
        return self.call_method("validate")

    # Do not allow setting of individual attributes
    def __setattr__(self, *args, **kwargs):
        raise Exception(
            "Please set all parameters at once via collision_detection.set_params()")

    # Override to call validate after parameter update
    def set_params(self, **kwargs):
        """
        Set the parameters for the collision detection

        See :ref:`Creating bonds when particles collide` for detailed instructions.


        Parameters
        ----------
        distance : :obj:`float`
            Distance below which a pair of particles is considered in the
            collision detection

        bond_type : :obj:`espressomd.interactions.BondedInteraction`
            Bond to add between the colliding particles

        vs_bond_type : :obj:`espressomd.interactions.BondedInteraction`
            Bond to add between virtual sites. If only one virtual site would be
            created, then the virtal site and the real particle, which does not hold
            this virtual site, will also connected via this bond.

        rate : :obj:`float`
            The rate at collision is accepted to form bonds.

        particle_type : array_like :obj:`int`
            The particle types, which will be checked for collision. At least two
            particle types need to be defined. If one want a collision of same type
            particles, the type can be added twice.

        vs_particle_type : aray_like :obj:`int`
            The virtual sites to be created between the colliding particles. This is
            be used to 'glue' the particles together, so that their relative orientation
            will remain the same. The entries of this array match the entries of the
            `particle_type` array, so that for each particle type there can be a
            different type of virtual site created.

        distance_vs_particle : aray_like :obj:`float`
            The relative distance (0<= d <= 1) at which the virtual site will be created.
            This distance needs to be set for each virtual site type created.

        bond_three_particles : :obj:`espressomd.interactions.BondedInteraction`
            First angular bond for the ``"bind_three_particles"`` mode. See
            user guide

        three_particle_binding_angle_resolution : :obj:`int`
            Resolution for the angular bonds (mode ``"bind_three_particles"``).
            Resolution+1 bonds are needed to accommodate the case of 180 degrees
            angles

        """

        # Convert bonds to bond ids
        for name in ["bond_type", "vs_bond_type"]:
            if name in kwargs:
                if isinstance(kwargs[name], BondedInteraction):
                    kwargs[name] = kwargs[name]._bond_id
        if not 'active' in kwargs and len(kwargs) > 0:
            kwargs['active'] = 1
        super().set_params(**kwargs)
        self.validate()
        handle_errors("Validation of collision detection failed")

    def get_parameter(self, name):
        """Gets a single parameter from the collision detection."""

        res = super().get_parameter(name)
        return self._convert_param(name, res)

    def get_params(self):
        """Returns the parameters of the collision detection as dict.

        """
        res = super().get_params()
        for k in res.keys():
            res[k] = self._convert_param(k, res[k])

        # Filter key-value pairs according to active mode
        return res

    def _convert_param(self, name, value):
        """
        Handles type conversion core -> python

        Bond types: int -> BondedInteraction
        mode: int -> string

        """
        # Py3: Cast from binary to normal string. Don't understand, why a
        # binary string can even occur, here, but it does.
        name = to_str(name)
        # Convert int mode parameter to string
        res = value

        # Convert bond parameters from bond ids to into BondedInteractions
        if name in ["bond_type", "vs_bond_type"]:
            if value == -1:  # Not defined
                res = None
            else:
                res = BondedInteractions()[value]
        return res

    # Pickle support
    def __reduce__(self):
        return _restore_collision_detection, (self.get_params(),)


def _restore_collision_detection(params):
    return CollisionDetection(**params)
