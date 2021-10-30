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
from libcpp cimport bool
include "myconfig.pxi"

import numpy as np
import collections

from .grid cimport box_geo
from . cimport integrate
from . import interactions
from . import integrate
from .actors import Actors
from . cimport cuda_init
from . import particle_data
from . import cuda_init
from . import code_info
from .utils cimport make_array_locked, make_Vector3d, Vector3d
from .thermostat import Thermostat
from .cellsystem import CellSystem
from .analyze import Analysis
from .galilei import GalileiTransform
from .constraints import Constraints
from .accumulators import AutoUpdateAccumulators
if LB_BOUNDARIES or LB_BOUNDARIES_GPU:
    from .lbboundaries import LBBoundaries
    from .ekboundaries import EKBoundaries
from .comfixed import ComFixed
from .utils cimport check_type_or_throw_except
from .utils import handle_errors, array_locked
IF VIRTUAL_SITES:
    from .virtual_sites import ActiveVirtualSitesHandle, VirtualSitesOff

IF COLLISION_DETECTION == 1:
    from .collision_detection import CollisionDetection


setable_properties = ["box_l", "min_global_cut", "periodicity", "time",
                      "time_step", "force_cap", "max_oif_objects"]
checkpointable_properties = ["max_oif_objects"]

if VIRTUAL_SITES:
    setable_properties.append("_active_virtual_sites_handle")
    checkpointable_properties.append("_active_virtual_sites_handle")


cdef bool _system_created = False

cdef class _Globals:
    def __getstate__(self):
        return {'box_l': self.box_l,
                'periodicity': self.periodicity,
                'min_global_cut': self.min_global_cut}

    def __setstate__(self, params):
        self.box_l = params['box_l']
        self.periodicity = params['periodicity']
        self.min_global_cut = params['min_global_cut']

    property box_l:
        """
        (3,) array_like of :obj:`float`:
            Dimensions of the simulation box

        """

        def __set__(self, _box_l):
            if len(_box_l) != 3:
                raise ValueError("Box length must be of length 3")
            mpi_set_box_length(make_Vector3d(_box_l))

        def __get__(self):
            return make_array_locked(< Vector3d > box_geo.length())

    property periodicity:
        """
        (3,) array_like of :obj:`bool`:
            System periodicity in ``[x, y, z]``, ``False`` for no periodicity
            in this direction, ``True`` for periodicity

        """

        def __set__(self, _periodic):
            if len(_periodic) != 3:
                raise ValueError(
                    f"periodicity must be of length 3, got length {len(_periodic)}")
            mpi_set_periodicity(_periodic[0], _periodic[1], _periodic[2])
            handle_errors("Error while assigning system periodicity")

        def __get__(self):
            periodicity = np.empty(3, dtype=type(True))
            for i in range(3):
                periodicity[i] = box_geo.periodic(i)
            return array_locked(periodicity)

    property min_global_cut:
        def __set__(self, _min_global_cut):
            mpi_set_min_global_cut(_min_global_cut)

        def __get__(self):
            global min_global_cut
            return min_global_cut

cdef class System:
    """The ESPResSo system class.

    .. note:: every attribute has to be declared at the class level.
              This means that methods cannot define an attribute by using
              ``self.new_attr = somevalue`` without declaring it inside this
              indentation level, either as method, property or reference.

    """

    cdef public:
        _globals
        part
        """:class:`espressomd.particle_data.ParticleList`"""
        non_bonded_inter
        """:class:`espressomd.interactions.NonBondedInteractions`"""
        bonded_inter
        """:class:`espressomd.interactions.BondedInteractions`"""
        cell_system
        """:class:`espressomd.cellsystem.CellSystem`"""
        thermostat
        """:class:`espressomd.thermostat.Thermostat`"""
        actors
        """:class:`espressomd.actors.Actors`"""
        analysis
        """:class:`espressomd.analyze.Analysis`"""
        galilei
        """:class:`espressomd.galilei.GalileiTransform`"""
        integrator
        """:class:`espressomd.integrate.IntegratorHandle`"""
        auto_update_accumulators
        """:class:`espressomd.accumulators.AutoUpdateAccumulators`"""
        constraints
        """:class:`espressomd.constraints.Constraints`"""
        lbboundaries
        """:class:`espressomd.lbboundaries.LBBoundaries`"""
        ekboundaries
        """:class:`espressomd.ekboundaries.EKBoundaries`"""
        collision_detection
        """:class:`espressomd.collision_detection.CollisionDetection`"""
        cuda_init_handle
        """:class:`espressomd.cuda_init.CudaInitHandle`"""
        comfixed
        """:class:`espressomd.comfixed.ComFixed`"""
        _active_virtual_sites_handle

    def __init__(self, **kwargs):
        global _system_created
        if not _system_created:
            if 'box_l' not in kwargs:
                raise ValueError("Required argument box_l not provided.")
            self._globals = _Globals()
            self.integrator = integrate.IntegratorHandle()
            System.__setattr__(self, "box_l", kwargs.pop("box_l"))
            for arg in kwargs:
                if arg in setable_properties:
                    System.__setattr__(self, arg, kwargs.get(arg))
                else:
                    raise ValueError(
                        f"Property {arg} can not be set via argument to System class.")
            self.actors = Actors()
            self.analysis = Analysis(self)
            self.auto_update_accumulators = AutoUpdateAccumulators()
            self.bonded_inter = interactions.BondedInteractions()
            self.cell_system = CellSystem()
            IF COLLISION_DETECTION == 1:
                self.collision_detection = CollisionDetection()
            self.comfixed = ComFixed()
            self.constraints = Constraints()
            IF CUDA:
                self.cuda_init_handle = cuda_init.CudaInitHandle()
            self.galilei = GalileiTransform()
            if LB_BOUNDARIES or LB_BOUNDARIES_GPU:
                self.lbboundaries = LBBoundaries()
                self.ekboundaries = EKBoundaries()
            self.non_bonded_inter = interactions.NonBondedInteractions()
            self.part = particle_data.ParticleList()
            self.thermostat = Thermostat()
            IF VIRTUAL_SITES:
                self._active_virtual_sites_handle = ActiveVirtualSitesHandle(
                    implementation=VirtualSitesOff())
            _system_created = True
        else:
            raise RuntimeError(
                "You can only have one instance of the system class at a time.")

    # __getstate__ and __setstate__ define the pickle interaction
    def __getstate__(self):
        odict = collections.OrderedDict()
        odict['_globals'] = System.__getattribute__(self, "_globals")
        odict['integrator'] = System.__getattribute__(self, "integrator")
        for property_ in checkpointable_properties:
            odict[property_] = System.__getattribute__(self, property_)
        odict['non_bonded_inter'] = System.__getattribute__(
            self, "non_bonded_inter")
        odict['bonded_inter'] = System.__getattribute__(self, "bonded_inter")
        odict['cell_system'] = System.__getattribute__(self, "cell_system")
        odict['part'] = System.__getattribute__(self, "part")
        odict['actors'] = System.__getattribute__(self, "actors")
        odict['analysis'] = System.__getattribute__(self, "analysis")
        odict['auto_update_accumulators'] = System.__getattribute__(
            self, "auto_update_accumulators")
        odict['comfixed'] = System.__getattribute__(self, "comfixed")
        odict['constraints'] = System.__getattribute__(self, "constraints")
        odict['galilei'] = System.__getattribute__(self, "galilei")
        IF LB_BOUNDARIES or LB_BOUNDARIES_GPU:
            odict['lbboundaries'] = System.__getattribute__(
                self, "lbboundaries")
        odict['thermostat'] = System.__getattribute__(self, "thermostat")
        IF COLLISION_DETECTION:
            odict['collision_detection'] = System.__getattribute__(
                self, "collision_detection")
        return odict

    def __setstate__(self, params):
        for property_ in params.keys():
            System.__setattr__(self, property_, params[property_])

    property box_l:
        """
        (3,) array_like of :obj:`float`:
            Dimensions of the simulation box

        """

        def __set__(self, _box_l):
            self._globals.box_l = _box_l

        def __get__(self):
            return self._globals.box_l

    property force_cap:
        """
        :obj:`float`:
            If > 0, the magnitude of the force on the particles
            are capped to this value.

        """

        def __get__(self):
            return self.integrator.force_cap

        def __set__(self, cap):
            self.integrator.force_cap = cap

    property periodicity:
        """
        (3,) array_like of :obj:`bool`:
            System periodicity in ``[x, y, z]``, ``False`` for no periodicity
            in this direction, ``True`` for periodicity

        """

        def __set__(self, _periodic):
            self._globals.periodicity = _periodic

        def __get__(self):
            return self._globals.periodicity

    property time:
        """
        Set the time in the simulation
        """

        def __set__(self, double sim_time):
            self.integrator.time = sim_time

        def __get__(self):
            return self.integrator.time

    property time_step:
        """
        Sets the time step for the integrator.
        """

        def __set__(self, double time_step):
            self.integrator.time_step = time_step

        def __get__(self):
            return self.integrator.time_step

    property max_cut_nonbonded:
        def __get__(self):
            return self.cell_system.max_cut_nonbonded

    property max_cut_bonded:
        def __get__(self):
            return self.cell_system.max_cut_bonded

    property min_global_cut:
        def __set__(self, _min_global_cut):
            self._globals.min_global_cut = _min_global_cut

        def __get__(self):
            return self._globals.min_global_cut

    IF VIRTUAL_SITES:
        property virtual_sites:
            def __set__(self, v):
                self._active_virtual_sites_handle.implementation = v

            def __get__(self):
                return self._active_virtual_sites_handle.implementation

    property max_oif_objects:
        """Maximum number of objects as per the object_in_fluid method.

        """

        def __get__(self):
            return max_oif_objects

        def __set__(self, v):
            mpi_set_max_oif_objects(v)

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

        if d_new < 0:
            raise ValueError("No negative lengths")
        if dir == "xyz":
            rescale_boxl(3, d_new)
        elif dir == "x" or dir == 0:
            rescale_boxl(0, d_new)
        elif dir == "y" or dir == 1:
            rescale_boxl(1, d_new)
        elif dir == "z" or dir == 2:
            rescale_boxl(2, d_new)
        else:
            raise ValueError(
                'Usage: change_volume_and_rescale_particles(<L_new>, [{ "x" | "y" | "z" | "xyz" }])')

    def volume(self):
        """Return box volume of the cuboid box.

        """

        return self.box_l[0] * self.box_l[1] * self.box_l[2]

    def distance(self, p1, p2):
        """Return the scalar distance between particles, between a particle
        and a point or between two points, respecting periodic boundaries.

        Parameters
        ----------
        p1 : :class:`~espressomd.particle_data.ParticleHandle` or (3,) array of :obj:`float`
            First particle or position.
        p2 : :class:`~espressomd.particle_data.ParticleHandle` or (3,) array of :obj:`float`
            Second particle or position.

        """
        res = self.distance_vec(p1, p2)
        return np.linalg.norm(res)

    def distance_vec(self, p1, p2):
        """Return the distance vector between particles, between a particle
        and a point or between two points, respecting periodic boundaries.

        Parameters
        ----------
        p1 : :class:`~espressomd.particle_data.ParticleHandle` or (3,) array of :obj:`float`
            First particle or position.
        p2 : :class:`~espressomd.particle_data.ParticleHandle` or (3,) array of :obj:`float`
            Second particle or position.

        """

        cdef Vector3d pos1
        if isinstance(p1, particle_data.ParticleHandle):
            pos1 = make_Vector3d(p1.pos)
        else:
            check_type_or_throw_except(
                p1, 3, float, "p1 must be a particle or 3 floats")
            pos1 = make_Vector3d(p1)
        cdef Vector3d pos2
        if isinstance(p2, particle_data.ParticleHandle):
            pos2 = make_Vector3d(p2.pos)
        else:
            check_type_or_throw_except(
                p2, 3, float, "p2 must be a particle or 3 floats")
            pos2 = make_Vector3d(p2)

        return make_array_locked(box_geo.get_mi_vector(pos2, pos1))

    def rotate_system(self, **kwargs):
        """Rotate the particles in the system about the center of mass.

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
        mpi_rotate_system(kwargs['phi'], kwargs['theta'], kwargs['alpha'])

    IF EXCLUSIONS:
        def auto_exclusions(self, distance):
            """Automatically adds exclusions between particles
            that are bonded.

            This only considers pair bonds.

            Parameters
            ----------
            distance : :obj:`int`
                Bond distance upto which the exclusions should be added.

            """
            auto_exclusions(distance)

    def setup_type_map(self, type_list=None):
        """
        For using ESPResSo conveniently for simulations in the grand canonical
        ensemble, or other purposes, when particles of certain types are created
        and deleted frequently. Particle ids can be stored in lists for each
        individual type and so random ids of particles of a certain type can be
        drawn. If you want ESPResSo to keep track of particle ids of a certain type
        you have to initialize the method by calling the setup function. After that
        ESPResSo will keep track of particle ids of that type.

        """
        if not hasattr(type_list, "__iter__"):
            raise ValueError("type_list has to be iterable.")

        for current_type in type_list:
            init_type_map(current_type)

    def number_of_particles(self, type=None):
        """
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

        """
        check_type_or_throw_except(type, 1, int, "type must be 1 int")
        number = number_of_particles_with_type(type)
        handle_errors("")
        return int(number)
