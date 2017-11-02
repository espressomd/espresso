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
if CONSTRAINTS == 1:
    from .constraints import Constraints

from .correlators import AutoUpdateCorrelators
from .observables import AutoUpdateObservables
if LB_BOUNDARIES or LB_BOUNDARIES_GPU:
    from .lbboundaries import LBBoundaries
from .ekboundaries import EKBoundaries
from .comfixed import ComFixed


import sys
import random  # for true random numbers from os.urandom()

setable_properties = ["box_l", "min_global_cut", "periodicity", "time",
                      "time_step", "timings", "force_cap"]

IF LEES_EDWARDS == 1:
    setable_properties.append("lees_edwards_offset")

cdef bool _system_created = False

cdef class System(object):
    """ The base class for espressomd.system.System().

    .. note:: every attribute has to be declared at the class level.
              This means that methods cannot define an attribute by using
              ``self.new_attr = somevalue`` without declaring it inside this
              indentation level, either as method, property or reference.

    """
    cdef public:
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
        auto_update_observables
        auto_update_correlators
        constraints
        lbboundaries
        ekboundaries
        __seed
        cuda_init_handle
        comfixed

    def __init__(self):
        global _system_created
        if (not _system_created):
            self.part = particle_data.ParticleList()
            self.non_bonded_inter = interactions.NonBondedInteractions()
            self.bonded_inter = interactions.BondedInteractions()
            self.cell_system = CellSystem()
            self.thermostat = Thermostat()
            self.minimize_energy = MinimizeEnergy()
            self.actors = Actors(_system=self)
            self.analysis = Analysis(self)
            self.galilei = GalileiTransform()
            self.integrator = integrate.Integrator()
            self.auto_update_observables = AutoUpdateObservables()
            self.auto_update_correlators = AutoUpdateCorrelators()
            if CONSTRAINTS:
                self.constraints = Constraints()
            if LB_BOUNDARIES or LB_BOUNDARIES_GPU:
                self.lbboundaries = LBBoundaries()
                self.ekboundaries = EKBoundaries()
            IF CUDA:
                self.cuda_init_handle = cuda_init.CudaInitHandle()

            self.comfixed = ComFixed()
            _system_created = True
        else:
            raise RuntimeError(
                "You can only have one instance of the system class at a time.")

    # __getstate__ and __setstate__ define the pickle interaction
    def __getstate__(self):
        odict = {}
        for property_ in setable_properties:
            odict[property_] = System.__getattribute__(self, property_)
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

            mpi_bcast_parameter(FIELD_BOXL)

        def __get__(self):
            return np.array([box_l[0], box_l[1], box_l[2]])

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
            return forcecap_get()

        def __set__(self, cap):
            forcecap_set(cap)

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
            periodicity = np.zeros(3)
            for i in range(3):
                if _periodic[i] != 1:
                    IF PARTIAL_PERIODIC:
                        pass
                    ELSE:
                        raise ValueError(
                            "The feature PARTIAL_PERIODIC needs to be activated in myconfig.hpp")
            for i in range(3):
                periodicity[i] = _periodic[i]
            periodic = 4 * _periodic[2] + 2 * _periodic[1] + _periodic[0]
            # first 3 bits of periodic determine the periodicity
            mpi_bcast_parameter(FIELD_PERIODIC)

        def __get__(self):
            periodicity = np.zeros(3)
            periodicity[0] = periodic % 2
            periodicity[1] = int(periodic / 2) % 2
            periodicity[2] = int(periodic / 4) % 2
            return periodicity

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

    property smaller_time_step:
        """
        Setting this property to a positive integer value turns on the multi-timestepping algorithm. The ratio :attr:`espressomd.system.System.time_step`/:attr:`espressomd.system.System.smaller_time_step` must be an integer.
        """
        def __set__(self, double _smaller_time_step):
            IF MULTI_TIMESTEP:
                global smaller_time_step
                if _smaller_time_step <= 0:
                    raise ValueError("Smaller time step must be positive")
                mpi_set_smaller_time_step(_smaller_time_step)

        def __get__(self):
            return smaller_time_step

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
            mpi_set_time_step(_time_step)

        def __get__(self):
            return time_step

    property timings:
        def __set__(self, int _timings):
            global timing_samples
            if _timings <= 0:
                timing_samples = 0
            else:
                timing_samples = _timings

        def __get__(self):
            return timing_samples

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
            global min_global_cut
            min_global_cut = _min_global_cut
            mpi_bcast_parameter(FIELD_MIN_GLOBAL_CUT)

        def __get__(self):
            return min_global_cut

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
            states[i] = " ".join(map(str, states_on_node_i))
        mpi_random_set_stat(states)

    property seed:
        """
        Sets the seed of the pseudo random number with a list of seeds which is as long as the number of used nodes.
        """
        
        def __set__(self, _seed):
            cdef vector[int] seed_array
            self.__seed = _seed
            if(isinstance(_seed, int) and n_nodes == 1):
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
                    states[i] = " ".join(
                        map(str, rng_state[i * _state_size_plus_one:(i + 1) * _state_size_plus_one]))
                mpi_random_set_stat(states)
            else:
                raise ValueError("Wrong # of args: Usage: 'random_number_generator_state \"<state(1)> ... <state(n_nodes*(state_size+1))>, where each <state(i)> is an integer. The state size of the PRNG can be obtained by calling _get_PRNG_state_size().")

        def __get__(self):
            rng_state = map(int, (mpi_random_get_stat().c_str()).split())
            return rng_state

    IF LEES_EDWARDS == 1:
        property lees_edwards_offset:
        # defines the lees edwards offset
            def __set__(self, double _lees_edwards_offset):

                if isinstance(_lees_edwards_offset, float):
                    global lees_edwards_offset
                    lees_edwards_offset = _lees_edwards_offset
                    #new_offset = _lees_edwards_offset
                    mpi_bcast_parameter(FIELD_LEES_EDWARDS_OFFSET)

                else:
                    raise ValueError("Wrong # of args! Usage: lees_edwards_offset { new_offset }")

            def __get__(self):
        # global lees_edwards_offset
                return lees_edwards_offset

    def change_volume_and_rescale_particles(self, d_new, dir="xyz"):
        """Change box size and rescale particle coordinates.

        Parameters:
        -----------
        d_new : float
                new box length
        dir : str, optional
                coordinate to work on, ``"x"``, ``"y"``, ``"z"`` or ``"xyz"`` for isotropic.
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
        """
        Return box volume of the cuboid box

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
