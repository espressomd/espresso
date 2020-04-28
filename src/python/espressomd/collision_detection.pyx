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
        # If no mode is specified at construction, use off.
        if "mode" not in kwargs:
            kwargs["mode"] = "off"
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
        mode : :obj:`str`, \{"off", "bind_centers", "bind_at_point_of_collision", "bind_three_particles", "glue_to_surface"\}
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

        bond_three_particles : :obj:`espressomd.interactions.BondedInteraction`
            First angular bond for the ``"bind_three_particles"`` mode. See
            user guide

        three_particle_binding_angle_resolution : :obj:`int`
            Resolution for the angular bonds (mode ``"bind_three_particles"``).
            Resolution+1 bonds are needed to accommodate the case of 180 degrees
            angles

        """

        if not ("mode" in kwargs):
            raise Exception(
                "Collision mode must be specified via the mode keyword argument")

        # Completeness of parameter set
        if not (set(kwargs.keys()) == set(
                self._params_for_mode(kwargs["mode"]))):
            raise Exception("Parameter set does not match mode. ",
                            kwargs["mode"], "requires ",
                            self._params_for_mode(kwargs["mode"]))

        # Mode
        kwargs["mode"] = self._int_mode[kwargs["mode"]]

        # Convert bonds to bond ids
        for name in ["bond_centers", "bond_vs", "bond_three_particle_binding"]:
            if name in kwargs:
                if isinstance(kwargs[name], BondedInteraction):
                    kwargs[name] = kwargs[name]._bond_id
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
        return {k: res[k] for k in self._params_for_mode(res["mode"])}

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
        if name == "mode":
            res = self._str_mode(value)

        # Convert bond parameters from bond ids to into BondedInteractions
        if name in ["bond_centers", "bond_vs", "bond_three_particle_binding"]:
            if value == -1:  # Not defined
                res = None
            else:
                res = BondedInteractions()[value]
        return res

    def _params_for_mode(self, mode):
        """The parameter names expected for a given collision mode

        """
        if mode == "off":
            return ("mode",)
        if mode == "bind_centers":
            return ("mode", "bond_centers", "distance")
        if mode == "bind_at_point_of_collision":
            return ("mode", "bond_centers", "bond_vs",
                    "part_type_vs", "distance", "vs_placement")
        if mode == "glue_to_surface":
            return ("mode", "bond_centers", "bond_vs", "part_type_vs",
                    "part_type_to_be_glued", "part_type_to_attach_vs_to",
                    "part_type_after_glueing", "distance",
                    "distance_glued_particle_to_vs")
        if mode == "bind_three_particles":
            return ("mode", "bond_centers", "distance", "bond_three_particles",
                    "three_particle_binding_angle_resolution")
        raise Exception("Mode not handled: " + mode.__str__())

    _int_mode = {
        "off": int(COLLISION_MODE_OFF),
        "bind_centers": int(COLLISION_MODE_BOND),
        "bind_at_point_of_collision": int(COLLISION_MODE_VS),
        "glue_to_surface": int(COLLISION_MODE_GLUE_TO_SURF),
        "bind_three_particles": int(COLLISION_MODE_BIND_THREE_PARTICLES)}

    def _str_mode(self, int_mode):
        """String mode name from int ones provided by the core

        """
        for key in self._int_mode:
            if self._int_mode[key] == int_mode:
                return key
        raise Exception("Unknown integer collision mode %d" % int_mode)

    # Pickle support
    def __reduce__(self):
        return _restore_collision_detection, (self.get_params(),)


def _restore_collision_detection(params):
    return CollisionDetection(**params)
