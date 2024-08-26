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

import numpy as np
import collections

# pylint: disable=unused-import
from . import accumulators
from . import analyze
from . import bond_breakage
from . import cell_system
from . import cuda_init
from . import collision_detection
from . import comfixed
from . import constraints
from . import electrostatics
from . import magnetostatics
from . import galilei
from . import interactions
from . import integrate
from . import lees_edwards
from . import particle_data
from . import thermostat
# pylint: enable=unused-import

from .code_features import has_features, assert_features
from .script_interface import script_interface_register, ScriptInterfaceHelper


@script_interface_register
class System(ScriptInterfaceHelper):
    """
    The ESPResSo system class.

    Attributes
    ----------
    analysis: :class:`espressomd.analyze.Analysis`
    auto_update_accumulators: :class:`espressomd.accumulators.AutoUpdateAccumulators`
    bond_breakage: :class:`espressomd.bond_breakage.BreakageSpecs`
    bonded_inter: :class:`espressomd.interactions.BondedInteractions`
    cell_system: :class:`espressomd.cell_system.CellSystem`
    collision_detection: :class:`espressomd.collision_detection.CollisionDetection`
    comfixed: :class:`espressomd.comfixed.ComFixed`
    constraints: :class:`espressomd.constraints.Constraints`
    cuda_init_handle: :class:`espressomd.cuda_init.CudaInitHandle`
    galilei: :class:`espressomd.galilei.GalileiTransform`
    integrator: :class:`espressomd.integrate.IntegratorHandle`
    lees_edwards: :class:`espressomd.lees_edwards.LeesEdwards`
    non_bonded_inter: :class:`espressomd.interactions.NonBondedInteractions`
    part: :class:`espressomd.particle_data.ParticleList`
    thermostat: :class:`espressomd.thermostat.Thermostat`
    box_l: (3,) array_like of :obj:`float`
        Dimensions of the simulation box.
    periodicity: (3,) array_like of :obj:`bool`
        System periodicity in ``[x, y, z]``, ``False`` for no periodicity
        in this direction, ``True`` for periodicity
    min_global_cut : :obj:`float`
        Minimal interaction cutoff.

    Methods
    -------
    setup_type_map()
        For using ESPResSo conveniently for simulations in the grand canonical
        ensemble, or other purposes, when particles of certain types are created
        and deleted frequently. Particle ids can be stored in lists for each
        individual type and so random ids of particles of a certain type can be
        drawn. If you want ESPResSo to keep track of particle ids of a certain type
        you have to initialize the method by calling the setup function. After that
        ESPResSo will keep track of particle ids of that type.

        Parameters
        ----------
        type_list : array_like of :obj:`int`
            Types to track.

    number_of_particles()
        Count the number of particles of a given type.

        Parameters
        ----------
        type : :obj:`int` (:attr:`~espressomd.particle_data.ParticleHandle.type`)
            Particle type to count the number for.

        Returns
        -------
        :obj:`int`
            The number of particles which have the given type.

        Raises
        ------
        RuntimeError
            If the particle ``type`` is not currently tracked by the system.
            To select which particle types are tracked, call :meth:`setup_type_map`.

    rotate_system()
        Rotate the particles in the system about the center of mass.

        If ``ROTATION`` is activated, the internal rotation degrees of
        freedom are rotated accordingly.

        Parameters
        ----------
        phi : :obj:`float`
            Angle between the z-axis and the rotation axis.
        theta : :obj:`float`
            Rotation of the axis around the y-axis.
        alpha : :obj:`float`
            How much to rotate

    """
    _so_name = "System::System"
    _so_creation_policy = "GLOBAL"
    _so_bind_methods = (
        "setup_type_map",
        "number_of_particles",
        "rotate_system")

    def __init__(self, **kwargs):
        if "sip" in kwargs:
            super().__init__(**kwargs)
            self._setup_atexit()
        else:
            super().__init__(_regular_constructor=True, **kwargs)
            if has_features("CUDA"):
                self.cuda_init_handle = cuda_init.CudaInitHandle()
            if has_features("WALBERLA"):
                self._lb = None
                self._ekcontainer = None

            # lock class
            self.call_method("lock_system_creation")
            self._setup_atexit()

        self._ase_interface = None

    def _setup_atexit(self):
        import atexit

        def session_shutdown():
            self.call_method("session_shutdown")
        atexit.register(session_shutdown)

    def __reduce__(self):
        so_callback, so_callback_args = super().__reduce__()
        return (System._restore_object,
                (so_callback, so_callback_args, self.__getstate__()))

    @classmethod
    def _restore_object(cls, so_callback, so_callback_args, state):
        so = so_callback(*so_callback_args)
        so.__setstate__(state)
        return so

    def __getstate__(self):
        checkpointable_properties = []
        if has_features("WALBERLA"):
            checkpointable_properties += ["_lb", "_ekcontainer"]

        odict = collections.OrderedDict()
        for property_name in checkpointable_properties:
            odict[property_name] = System.__getattribute__(self, property_name)
        if self._ase_interface is not None:
            odict["_ase_interface"] = self._ase_interface.__getstate__()
        return odict

    def __setstate__(self, params):
        # initialize Python-only members
        if "_ase_interface" in params:
            from espressomd.plugins.ase import ASEInterface
            self.ase = ASEInterface(**params.pop("_ase_interface"))
        for property_name in params.keys():
            System.__setattr__(self, property_name, params[property_name])
        # note: several members can only be instantiated once
        if has_features("WALBERLA"):
            if self._lb is not None:
                lb, self._lb = self._lb, None
                self.lb = lb
            if self._ekcontainer is not None:
                ekcontainer, self._ekcontainer = self._ekcontainer, None
                self.ekcontainer = ekcontainer
        self.call_method("lock_system_creation")

    @property
    def force_cap(self):
        """
        If > 0, the magnitude of the force on the particles
        are capped to this value.

        Type: :obj:`float`

        """
        return self.integrator.force_cap

    @force_cap.setter
    def force_cap(self, value):
        self.integrator.force_cap = value

    @property
    def time(self):
        """
        Total simulation time.

        Type: :obj:`float`

        """
        return self.integrator.time

    @time.setter
    def time(self, value):
        self.integrator.time = value

    @property
    def time_step(self):
        """
        MD time step.

        Type: :obj:`float`

        """
        return self.integrator.time_step

    @time_step.setter
    def time_step(self, value):
        self.integrator.time_step = value

    @property
    def max_cut_nonbonded(self):
        """
        Maximal cutoff for non-bonded interactions.

        Type: :obj:`float`

        """
        return self.cell_system.max_cut_nonbonded

    @property
    def max_cut_bonded(self):
        """
        Maximal cutoff for bonded interactions.

        Type: :obj:`float`

        """
        return self.cell_system.max_cut_bonded

    @property
    def lb(self):
        """
        LB solver.

        """
        assert_features("WALBERLA")
        return self._lb

    @lb.setter
    def lb(self, lb):
        assert_features("WALBERLA")
        if lb != self._lb:
            if self._lb is not None:
                self._lb.call_method("deactivate")
                self._lb = None
            if lb is not None:
                lb.call_method("activate")
                self._lb = lb

    @property
    def ase(self):
        return self._ase_interface

    @ase.setter
    def ase(self, ase):
        ase.register_system(self)
        self._ase_interface = ase

    @property
    def ekcontainer(self):
        """
        EK system (diffusion-advection-reaction models).

        Type: :class:`espressomd.electrokinetics.EKContainer`

        """
        assert_features("WALBERLA")
        return self._ekcontainer

    @ekcontainer.setter
    def ekcontainer(self, ekcontainer):
        assert_features("WALBERLA")
        if ekcontainer != self._ekcontainer:
            if self._ekcontainer is not None:
                self._ekcontainer.call_method("deactivate")
                self._ekcontainer = None
            if ekcontainer is not None:
                ekcontainer.call_method("activate")
                self._ekcontainer = ekcontainer

    def change_volume_and_rescale_particles(self, d_new, dir="xyz"):
        """Change box size and rescale particle coordinates.

        Parameters
        ----------
        d_new : :obj:`float`
            New box length
        dir : :obj:`str`, optional
            Coordinate to work on, ``"x"``, ``"y"``, ``"z"`` or ``"xyz"`` for isotropic.
            Isotropic assumes a cubic box.

        """

        coord = {"x": 0, "y": 1, "z": 2, 0: 0, 1: 1, 2: 2, "xyz": 3}.get(dir)
        if coord is None:
            raise ValueError(
                'Usage: change_volume_and_rescale_particles(<L_new>, [{ "x" | "y" | "z" | "xyz" }])')
        self.call_method("rescale_boxl", length=d_new, coord=coord)

    def volume(self):
        """Return volume of the cuboid box.

        """

        return float(np.prod(self.box_l))

    def distance(self, p1, p2):
        """Return the scalar distance between particles, between a particle
        and a point or between two points, respecting periodic boundaries.

        Parameters
        ----------
        p1 : :class:`~espressomd.particle_data.ParticleHandle` or (3,) array_like of :obj:`float`
            First particle or position.
        p2 : :class:`~espressomd.particle_data.ParticleHandle` or (3,) array_like of :obj:`float`
            Second particle or position.

        """
        return np.linalg.norm(self.distance_vec(p1, p2))

    def distance_vec(self, p1, p2):
        """Return the distance vector between particles, between a particle
        and a point or between two points, respecting periodic boundaries.

        Parameters
        ----------
        p1 : :class:`~espressomd.particle_data.ParticleHandle` or (3,) array_like of :obj:`float`
            First particle or position.
        p2 : :class:`~espressomd.particle_data.ParticleHandle` or (3,) array_like of :obj:`float`
            Second particle or position.

        """

        if isinstance(p1, particle_data.ParticleHandle):
            pos1 = p1.pos_folded
        else:
            pos1 = p1
        if isinstance(p2, particle_data.ParticleHandle):
            pos2 = p2.pos_folded
        else:
            pos2 = p2

        return self.call_method("distance_vec", pos1=pos1, pos2=pos2)

    def velocity_difference(self, p1, p2):
        """
        Return the velocity difference between two particles,
        considering Lees-Edwards boundary conditions, if active.

        Parameters
        ----------
        p1 : :class:`~espressomd.particle_data.ParticleHandle`
        p2 : :class:`~espressomd.particle_data.ParticleHandle`

        """

        return self.call_method("velocity_difference", pos1=p1.pos_folded,
                                pos2=p2.pos_folded, v1=p1.v, v2=p2.v)

    def auto_exclusions(self, distance):
        """
        Add exclusions between particles that are bonded.

        This only considers pair bonds.

        Requires feature ``EXCLUSIONS``.

        Parameters
        ----------
        distance : :obj:`int`
            Bond distance up to which the exclusions should be added.

        """
        assert_features("EXCLUSIONS")
        self.part.auto_exclusions(distance=distance)
