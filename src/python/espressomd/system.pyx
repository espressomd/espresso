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
from __future__ import print_function, absolute_import
from libcpp cimport bool
include "myconfig.pxi"

from globals cimport *
import numpy as np
import collections

from . cimport integrate
from . import interactions
from . import integrate
from .actors import Actors
from . cimport cuda_init
from . import particle_data
from . import cuda_init
from . import code_info
from .utils cimport numeric_limits
from .thermostat import Thermostat
from .cellsystem import CellSystem
from .minimize_energy import MinimizeEnergy
from .analyze import Analysis
from .galilei import GalileiTransform
from .constraints import Constraints

from .accumulators import AutoUpdateAccumulators
if LB_BOUNDARIES or LB_BOUNDARIES_GPU:
    from .lbboundaries import LBBoundaries
from .ekboundaries import EKBoundaries
from .comfixed import ComFixed
from globals cimport max_seen_particle
from .globals import Globals
from espressomd.utils import array_locked, is_valid_type
from espressomd.virtual_sites import ActiveVirtualSitesHandle, VirtualSitesOff

IF COLLISION_DETECTION == 1:
    from .collision_detection import CollisionDetection

import sys
import random  # for true random numbers from os.urandom()
cimport tuning


setable_properties = ["box_l", "min_global_cut", "periodicity", "time",
                      "time_step", "timings", "force_cap"]

if VIRTUAL_SITES:
    setable_properties.append("_active_virtual_sites_handle")


if OIF_GLOBAL_FORCES:
    setable_properties.append("max_oif_objects")

cdef bool _system_created = False

cdef class System(object):
    """ The base class for espressomd.system.System().

    .. note:: every attribute has to be declared at the class level.
              This means that methods cannot define an attribute by using
              ``self.new_attr = somevalue`` without declaring it inside this
              indentation level, either as method, property or reference.

    """
    cdef public:
        globals
        part
        non_bonded_inter
        bonded_inter
        cell_system
        thermostat
        minimize_energy
        actors
        analysis
        galilei
        integrator
        auto_update_accumulators
        constraints
        lbboundaries
        ekboundaries
        collision_detection
        __seed
        cuda_init_handle
        comfixed
        _active_virtual_sites_handle

    def __init__(self, **kwargs):
        global _system_created
        if (not _system_created):
            self.globals = Globals()
            if 'box_l' not in kwargs:
                raise ValueError("Required argument box_l not provided.")
            System.__setattr__(self, "box_l", kwargs.get("box_l"))
            del kwargs["box_l"]
            for arg in kwargs:
                if arg in setable_properties:
                    System.__setattr__(self, arg, kwargs.get(arg))
                else:
                    raise ValueError("Property {} can not be set via argument to System class.".format(arg))
            self.actors = Actors()
            self.analysis = Analysis(self)
            self.auto_update_accumulators = AutoUpdateAccumulators()
            self.bonded_inter = interactions.BondedInteractions()
            self.cell_system = CellSystem()
            IF COLLISION_DETECTION==1:
                self.collision_detection = CollisionDetection()
            self.comfixed = ComFixed()
            self.constraints = Constraints()
            IF CUDA:
                self.cuda_init_handle = cuda_init.CudaInitHandle()
            self.galilei = GalileiTransform()
            self.integrator = integrate.Integrator()
            if LB_BOUNDARIES or LB_BOUNDARIES_GPU:
                self.lbboundaries = LBBoundaries()
                self.ekboundaries = EKBoundaries()
            self.minimize_energy = MinimizeEnergy()
            self.non_bonded_inter = interactions.NonBondedInteractions()
            self.part = particle_data.ParticleList()
            self.thermostat = Thermostat()
            IF VIRTUAL_SITES:
                self._active_virtual_sites_handle=ActiveVirtualSitesHandle(implementation=VirtualSitesOff())
            _system_created = True
        else:
            raise RuntimeError(
                "You can only have one instance of the system class at a time.")

    # __getstate__ and __setstate__ define the pickle interaction
    def __getstate__(self):
        odict = collections.OrderedDict()
        odict['globals'] = System.__getattribute__(self, "globals")
        for property_ in setable_properties:
            if not hasattr(self.globals, property_):
                odict[property_] = System.__getattribute__(self, property_)
        odict['non_bonded_inter'] = System.__getattribute__(self, "non_bonded_inter")
        odict['bonded_inter'] = System.__getattribute__(self, "bonded_inter")
        odict['part'] = System.__getattribute__(self, "part")
        odict['actors'] = System.__getattribute__(self, "actors")
        odict['analysis'] = System.__getattribute__(self, "analysis")
        odict['auto_update_accumulators'] = System.__getattribute__(self, "auto_update_accumulators")
        odict['cell_system'] = System.__getattribute__(self, "cell_system")
        odict['comfixed'] = System.__getattribute__(self, "comfixed")
        odict['constraints'] = System.__getattribute__(self, "constraints")
        odict['galilei'] = System.__getattribute__(self, "galilei")
        odict['integrator'] = System.__getattribute__(self, "integrator")
        IF LB_BOUNDARIES or LB_BOUNDARIES_GPU:
            odict['lbboundaries'] = System.__getattribute__(self, "lbboundaries")
        odict['minimize_energy'] = System.__getattribute__(self, "minimize_energy")
        odict['thermostat'] = System.__getattribute__(self, "thermostat")
        return odict

    def __setstate__(self, params):
        for property_ in params.keys():
            System.__setattr__(self, property_, params[property_])

    property box_l:
        """
        Array like, list of three floats

        """

        def __set__(self, _box_l):
            if len(_box_l) != 3:
                raise ValueError("Box length must be of length 3")
            for i in range(3):
                if _box_l[i] <= 0:
                    raise ValueError(
                        "Box length must be > 0 in all directions")
                box_l[i] = _box_l[i]

            self.globals.box_l = box_l

        def __get__(self):
            return self.globals.box_l

    property integ_switch:
        def __get__(self):
            return integ_switch

    property force_cap:
        """
        If > 0, the magnitude of the force on the particles
        are capped to this value.

        type : float

        """
        def __get__(self):
            return self.globals.force_cap

        def __set__(self, cap):
            self.globals.force_cap = cap

    property periodicity:
        """
        list of three integers
        [x, y, z]
        zero for no periodicity in this direction
        one for periodicity

        """

        def __set__(self, _periodic):
            global periodic
            if len(_periodic) != 3:
                raise ValueError(
                    "periodicity must be of length 3, got length " + str(len(_periodic)))
            for i in range(3):
                if _periodic[i] != 1:
                    IF PARTIAL_PERIODIC:
                        pass
                    ELSE:
                        raise ValueError(
                            "The feature PARTIAL_PERIODIC needs to be activated in myconfig.hpp")
            self.globals.periodicity = _periodic

        def __get__(self):
            return self.globals.periodicity

    property time:
        """
        Set the time in the simulation 
        """
        def __set__(self, double _time):
            if _time < 0:
                raise ValueError("Simulation time must be >= 0")
            global sim_time
            sim_time = _time
            mpi_bcast_parameter(FIELD_SIMTIME)

        def __get__(self):
            global sim_time
            return sim_time

    property time_step:
        """
        Sets the time step for the integrator. 
        """
        def __set__(self, double _time_step):
            IF LB:
                global lbpar
            IF LB_GPU:
                global lbpar_gpu
            if _time_step <= 0:
                raise ValueError("Time Step must be positive")
            IF LB:
                if lbpar.tau >= 0.0 and _time_step < lbpar.tau:
                    raise ValueError(
                        "Time Step (" + str(time_step) + ") must be > LB_time_step (" + str(lbpar.tau) + ")")
            IF LB_GPU:
                if (lbpar_gpu.tau >= 0.0 and
                        lbpar_gpu.tau - _time_step > numeric_limits[float].epsilon() * abs(lbpar_gpu.tau + _time_step)):
                    raise ValueError(
                        "Time Step (" + str(time_step) + ") must be > LB_time_step (" + str(lbpar_gpu.tau) + ")")
            self.globals.time_step = _time_step

        def __get__(self):
            return self.globals.time_step

    property timings:
        def __set__(self, int _timings):
            self.globals.timings = _timings

        def __get__(self):
            return self.globals.timings

    property max_cut_nonbonded:
        def __get__(self):
            return max_cut_nonbonded

    property lattice_switch:
        def __get__(self):
            return lattice_switch

    property max_cut_bonded:
        def __get__(self):
            return max_cut_bonded

    property min_global_cut:
        def __set__(self, _min_global_cut):
            self.globals.min_global_cut = _min_global_cut

        def __get__(self):
            return self.globals.min_global_cut

    def _get_PRNG_state_size(self):
        """
        Returns the state of the pseudo random number generator.
        """
        
        return get_state_size_of_generator()

    def set_random_state_PRNG(self):
        """
        Sets the state of the pseudo random number generator using real random numbers.
        """
        
        _state_size_plus_one = self._get_PRNG_state_size() + 1
        states = string_vec(n_nodes)
        rng = random.SystemRandom()  # true RNG that uses os.urandom()
        for i in range(n_nodes):
            states_on_node_i = []
            for j in range(_state_size_plus_one + 1):
                states_on_node_i.append(
                    rng.randint(0, numeric_limits[int].max()))
            states[i] = (" ".join(map(str, states_on_node_i))).encode('utf-8')
        mpi_random_set_stat(states)

    property seed:
        """
        Sets the seed of the pseudo random number with a list of seeds which is as long as the number of used nodes.
        """

        def __set__(self, _seed):
            cdef vector[int] seed_array
            self.__seed = _seed
            if(is_valid_type(_seed, int) and n_nodes == 1):
                seed_array.resize(1)
                seed_array[0] = int(_seed)
                mpi_random_seed(0, seed_array)
            elif(hasattr(_seed, "__iter__")):
                if(len(_seed) < n_nodes or len(_seed) > n_nodes):
                    raise ValueError(
                        "The list needs to contain one seed value per node")
                seed_array.resize(len(_seed))
                for i in range(len(_seed)):
                    seed_array[i] = int(_seed[i])
                mpi_random_seed(n_nodes, seed_array)
            else:
                raise ValueError(
                    "The seed has to be an integer or a list of integers with one integer per node")

        def __get__(self):
            return self.__seed

    property random_number_generator_state:
        """Sets the random number generator state in the core. this is of interest for deterministic checkpointing
        """
        
        def __set__(self, rng_state):
            _state_size_plus_one = self._get_PRNG_state_size() + 1
            if(len(rng_state) == n_nodes * _state_size_plus_one):
                states = string_vec(n_nodes)
                for i in range(n_nodes):
                    states[i] = (" ".join(map(str,
                                 rng_state[i*_state_size_plus_one:(i+1)*_state_size_plus_one])
                                 )).encode('utf-8')
                mpi_random_set_stat(states)
            else:
                raise ValueError("Wrong # of args: Usage: 'random_number_generator_state \"<state(1)> ... <state(n_nodes*(state_size+1))>, where each <state(i)> is an integer. The state size of the PRNG can be obtained by calling _get_PRNG_state_size().")

        def __get__(self):
            rng_state = list(map(int, (mpi_random_get_stat().c_str()).split()))
            return rng_state

    IF VIRTUAL_SITES:
        property virtual_sites:
            def __set__(self,v):
                self._active_virtual_sites_handle.implementation=v
            def __get__(self):
                return self._active_virtual_sites_handle.implementation

    
    IF OIF_GLOBAL_FORCES:
        property max_oif_objects:
            """Maximum number of objects as per the object_in_fluid method.

            """
            
            def __get__(self):
                return max_oif_objects

            def __set__(self,v):
                global max_oif_objects
                max_oif_objects=v
                mpi_bcast_parameter(FIELD_MAX_OIF_OBJECTS)

    
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
        """Return the scalar distance between the particles, respecting periodic boundaries.

        """
        res = self.distance_vec(p1, p2)
        return np.sqrt(res[0]**2 + res[1]**2 + res[2]**2)

    def distance_vec(self, p1, p2):
        """Return the distance vector between the particles, respecting periodic boundaries.

        """
        cdef double[3] res, a, b
        a = p1.pos
        b = p2.pos

        get_mi_vector(res, b, a)
        return np.array((res[0], res[1], res[2]))

    def rotate_system(self, **kwargs):
        """Rotate the particles in the system about the center of mass.

           If ROTATION is activated, the internal rotation degrees of
           freedom are rotated accordingly.

        Parameters
        ----------
        phi : :obj:`float`
                Angle between the z-axis and the roation axis.
        theta : :obj:`float`
                Rotaton of the axis around the y-axis.
        alpha : :obj:`float`
                How much to rotate

        """
        rotate_system(kwargs['phi'], kwargs['theta'], kwargs['alpha'])

    IF EXCLUSIONS:
        def auto_exclusions(self, distance):
            """Automatically adds exclusions between particles
            that are bonded.

            This only considers pair bonds.

            Parameters
            ----------
            distance : :obj:`int`
                       Bond distance upto which the exlucsions should be added.

            """
            auto_exclusions(distance)


    def _is_valid_type(self, current_type):
        return (not (isinstance(current_type, int) or current_type < 0 or current_type > globals.max_seen_particle_type))


    def check_valid_type(self, current_type):
        if self._is_valid_type(current_type):
            raise ValueError("type", current_type, "does not exist!")


    def setup_type_map(self, type_list=None):
        """
        For using Espresso conveniently for simulations in the grand canonical
        ensemble, or other purposes, when particles of certain types are created
        and deleted frequently. Particle ids can be stored in lists for each
        individual type and so random ids of particles of a certain type can be
        drawn. If you want Espresso to keep track of particle ids of a certain type
        you have to initialize the method by calling the setup function. After that
        Espresso will keep track of particle ids of that type.

        """
        if not hasattr(type_list, "__iter__"):
            raise ValueError("type_list has to be iterable.")

        for current_type in type_list:
            init_type_map(current_type)

    def number_of_particles(self, type=None):
        """
        Parameters
        ----------
        current_type : :obj:`int` (:attr:`espressomd.particle_data.ParticleHandle.type`)
                       Particle type to count the number for. 

        Returns
        -------
        :obj:`int`
            The number of particles which have the given type.

        """
        self.check_valid_type( type)
        number=number_of_particles_with_type(type)
        return int(number)

    def find_particle(self, type=None):
        """
        The command will return a randomly chosen particle id, for a particle of
        the given type.
        
        """
        self.check_valid_type(type)
        pid=get_random_p_id(type)
        return int(pid)
