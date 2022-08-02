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
include "myconfig.pxi"

import numpy as np
import collections

from . import accumulators
from . import actors
from . import analyze
from . import bond_breakage
from . import cell_system
from . import cuda_init
from . import collision_detection
from . import comfixed
from . import constraints
from . import galilei
from . import interactions
from . import integrate
from . import lb
from . import EKSpecies
from . import lees_edwards
from . import particle_data
from . import thermostat
from . import virtual_sites

from .__init__ import has_features, assert_features
from .grid cimport box_geo
from .utils cimport Vector3d
from . cimport utils
from . import utils


_system_created = False


cdef class _BoxGeometry:
    """
    Wrapper class required for technical reasons only.

    When reloading from a checkpoint file, the box length, periodicity, and
    global cutoff must be set before anything else. Due to how pickling works,
    this can only be achieved by encapsulating them in a member object of the
    System class, and adding that object as the first element of the ordered
    dict that is used during serialization. When the System class is reloaded,
    the ordered dict is walked through and objects are deserialized in the same
    order. Since many objects depend on the box length, the `_BoxGeometry` has
    to be deserialized first. This guarantees the box geometry is already set
    in the core before e.g. particles and bonds are deserialized.

    """

    def __getstate__(self):
        return {'box_l': self.box_l,
                'periodicity': self.periodicity,
                'min_global_cut': self.min_global_cut}

    def __setstate__(self, params):
        self.box_l = params['box_l']
        self.periodicity = params['periodicity']
        self.min_global_cut = params['min_global_cut']

    property box_l:
        def __set__(self, box_l):
            utils.check_type_or_throw_except(
                box_l, 3, float, "box_l must be an array_like of 3 floats")
            mpi_set_box_length(utils.make_Vector3d(box_l))
            utils.handle_errors("Exception while updating the box length")

        def __get__(self):
            return utils.make_array_locked(box_geo.length())

    property periodicity:
        def __set__(self, periodic):
            utils.check_type_or_throw_except(
                periodic, 3, type(True), "periodicity must be an array_like of 3 bools")
            mpi_set_periodicity(periodic[0], periodic[1], periodic[2])
            utils.handle_errors("Exception while assigning system periodicity")

        def __get__(self):
            periodicity = np.empty(3, dtype=type(True))
            for i in range(3):
                periodicity[i] = box_geo.periodic(i)
            return utils.array_locked(periodicity)

    property min_global_cut:
        def __set__(self, min_global_cut):
            mpi_set_min_global_cut(min_global_cut)

        def __get__(self):
            return get_min_global_cut()


cdef class System:
    """The ESPResSo system class.

    .. note:: every attribute has to be declared at the class level.
              This means that methods cannot define an attribute by using
              ``self.new_attr = somevalue`` without declaring it inside this
              indentation level, either as method, property or reference.

    """

    cdef public:
        _box_geo
        part
        """:class:`espressomd.particle_data.ParticleList`"""
        non_bonded_inter
        """:class:`espressomd.interactions.NonBondedInteractions`"""
        bonded_inter
        """:class:`espressomd.interactions.BondedInteractions`"""
        cell_system
        """:class:`espressomd.cell_system.CellSystem`"""
        thermostat
        """:class:`espressomd.thermostat.Thermostat`"""
        actors
        """:class:`espressomd.actors.Actors`"""
        analysis
        """:class:`espressomd.analyze.Analysis`"""
        bond_breakage
        """:class:`espressomd.bond_breakage.BreakageSpecs`"""
        galilei
        """:class:`espressomd.galilei.GalileiTransform`"""
        integrator
        """:class:`espressomd.integrate.IntegratorHandle`"""
        auto_update_accumulators
        """:class:`espressomd.accumulators.AutoUpdateAccumulators`"""
        constraints
        """:class:`espressomd.constraints.Constraints`"""
        ekcontainer
        """:class:`espressomd.EKSpecies.EKContainer`"""
        ekreactions
        """:class:`espressomd.EKSpecies.EKReactions`"""
        lees_edwards
        """:class:`espressomd.lees_edwards.LeesEdwards`"""
        collision_detection
        """:class:`espressomd.collision_detection.CollisionDetection`"""
        cuda_init_handle
        """:class:`espressomd.cuda_init.CudaInitHandle`"""
        comfixed
        """:class:`espressomd.comfixed.ComFixed`"""
        _active_virtual_sites_handle

    def __init__(self, **kwargs):
        if _system_created:
            raise RuntimeError(
                "You can only have one instance of the system class at a time.")
        if 'box_l' not in kwargs:
            raise ValueError("Required argument 'box_l' not provided.")

        setable_properties = ["box_l", "min_global_cut", "periodicity", "time",
                              "time_step", "force_cap", "max_oif_objects"]
        if has_features("VIRTUAL_SITES"):
            setable_properties.append("_active_virtual_sites_handle")

        self._box_geo = _BoxGeometry()
        self.integrator = integrate.IntegratorHandle()
        System.__setattr__(self, "box_l", kwargs.pop("box_l"))
        for arg in kwargs:
            if arg not in setable_properties:
                raise ValueError(
                    f"Property '{arg}' can not be set via argument to System class.")
            System.__setattr__(self, arg, kwargs.get(arg))
        self.actors = actors.Actors()
        self.analysis = analyze.Analysis(self)
        self.auto_update_accumulators = accumulators.AutoUpdateAccumulators()
        self.bonded_inter = interactions.BondedInteractions()
        self.cell_system = cell_system.CellSystem()
        self.bond_breakage = bond_breakage.BreakageSpecs()
        if has_features("COLLISION_DETECTION"):
            self.collision_detection = collision_detection.CollisionDetection(
                mode="off")
        self.comfixed = comfixed.ComFixed()
        self.constraints = constraints.Constraints()
        if has_features("CUDA"):
            self.cuda_init_handle = cuda_init.CudaInitHandle()
        if has_features("LB_WALBERLA"):
            self.ekcontainer = EKSpecies.EKContainer()
            self.ekreactions = EKSpecies.EKReactions()
        self.galilei = galilei.GalileiTransform()
        self.lees_edwards = lees_edwards.LeesEdwards()
        self.non_bonded_inter = interactions.NonBondedInteractions()
        self.part = particle_data.ParticleList()
        self.thermostat = thermostat.Thermostat()
        if has_features("VIRTUAL_SITES"):
            self._active_virtual_sites_handle = virtual_sites.ActiveVirtualSitesHandle(
                implementation=virtual_sites.VirtualSitesOff())

        # lock class
        global _system_created
        _system_created = True

    # __getstate__ and __setstate__ define the pickle interaction
    def __getstate__(self):
        checkpointable_properties = ["_box_geo", "integrator"]
        if has_features("VIRTUAL_SITES"):
            checkpointable_properties.append("_active_virtual_sites_handle")
        checkpointable_properties += [
            "non_bonded_inter", "bonded_inter", "cell_system", "lees_edwards",
            "part", "analysis", "auto_update_accumulators", "comfixed",
            "constraints", "galilei", "bond_breakage", "max_oif_objects"
        ]
        if has_features("COLLISION_DETECTION"):
            checkpointable_properties.append("collision_detection")
        checkpointable_properties += ["actors", "thermostat"]
        if has_features("LB_WALBERLA"):
            checkpointable_properties += ["ekcontainer", "ekreactions"]

        odict = collections.OrderedDict()
        for property_name in checkpointable_properties:
            odict[property_name] = System.__getattribute__(self, property_name)
        if has_features("LB_WALBERLA"):
            odict["_lb_vtk_registry"] = lb._vtk_registry
            # TODO walberla
            # odict["_ek_vtk_registry"] = EKSpecies._ek_vtk_registry
        return odict

    def __setstate__(self, params):
        vtk_registry = None
        if has_features("LB_WALBERLA"):
            lb_vtk_registry = params.pop("_lb_vtk_registry")
            # TODO walberla
            # ek_vtk_registry = params.pop("_ek_vtk_registry")
        for property_name in params.keys():
            System.__setattr__(self, property_name, params[property_name])
        if has_features("LB_WALBERLA"):
            lb._vtk_registry = lb_vtk_registry
            # TODO walberla
            # EKSpecies._vtk_registry = ek_vtk_registry

    property box_l:
        """
        (3,) array_like of :obj:`float`:
            Dimensions of the simulation box

        """

        def __set__(self, value):
            self._box_geo.box_l = value

        def __get__(self):
            return self._box_geo.box_l

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

        def __set__(self, value):
            self._box_geo.periodicity = value

        def __get__(self):
            return self._box_geo.periodicity

    property time:
        """
        Set the time in the simulation.

        """

        def __set__(self, double sim_time):
            self.integrator.time = sim_time

        def __get__(self):
            return self.integrator.time

    property time_step:
        """
        Set the time step for the integrator.

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
        def __set__(self, value):
            self._box_geo.min_global_cut = value

        def __get__(self):
            return self._box_geo.min_global_cut

    property virtual_sites:
        """
        Set the virtual site implementation.

        Requires feature ``VIRTUAL_SITES``.

        """

        def __set__(self, v):
            assert_features("VIRTUAL_SITES")
            self._active_virtual_sites_handle.implementation = v

        def __get__(self):
            assert_features("VIRTUAL_SITES")
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
            pos1 = utils.make_Vector3d(p1.pos_folded)
        else:
            utils.check_type_or_throw_except(
                p1, 3, float, "p1 must be a particle or 3 floats")
            pos1 = utils.make_Vector3d(p1)
        cdef Vector3d pos2
        if isinstance(p2, particle_data.ParticleHandle):
            pos2 = utils.make_Vector3d(p2.pos_folded)
        else:
            utils.check_type_or_throw_except(
                p2, 3, float, "p2 must be a particle or 3 floats")
            pos2 = utils.make_Vector3d(p2)

        return utils.make_array_locked(box_geo.get_mi_vector(pos2, pos1))

    def velocity_difference(self, p1, p2):
        """
        Return the velocity difference between two particles,
        considering Lees-Edwards boundary conditions, if active.

        Parameters
        ----------
        p1 : :class:`~espressomd.particle_data.ParticleHandle`
        p2 : :class:`~espressomd.particle_data.ParticleHandle`

        """

        cdef Vector3d pos1 = utils.make_Vector3d(p1.pos_folded)
        cdef Vector3d pos2 = utils.make_Vector3d(p2.pos_folded)

        cdef Vector3d v1 = utils.make_Vector3d(p1.v)
        cdef Vector3d v2 = utils.make_Vector3d(p2.v)
        cdef Vector3d vd = box_geo.velocity_difference(pos2, pos1, v2, v1)

        return utils.make_array_locked(vd)

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

    def auto_exclusions(self, distance):
        """Automatically adds exclusions between particles
        that are bonded.

        This only considers pair bonds.

        Requires feature ``EXCLUSIONS``.

        Parameters
        ----------
        distance : :obj:`int`
            Bond distance upto which the exclusions should be added.

        """
        IF EXCLUSIONS:
            auto_exclusions(distance)
        ELSE:
            assert_features("EXCLUSIONS")

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
        utils.check_type_or_throw_except(type, 1, int, "type must be 1 int")
        return number_of_particles_with_type(type)
