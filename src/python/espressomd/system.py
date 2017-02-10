from __future__ import print_function
from espressomd.highlander import highlander
from espressomd import _system


import numpy as np

from . import integrate
from . import interactions
from . import integrate
from .actors import Actors
from . import cuda_init
from . import particle_data
from . import cuda_init
from . import code_info
from .thermostat import Thermostat
from .cellsystem import CellSystem
from .minimize_energy import MinimizeEnergy
from .analyze import Analysis
from .galilei import GalileiTransform
if "CONSTRAINTS" in code_info.features():
    from .constraints import Constraints

from .correlators import AutoUpdateCorrelators
from .observables import AutoUpdateObservables
from .lbboundaries import LBBoundaries
from .ekboundaries import EKBoundaries

import sys
import random  # for true random numbers from os.urandom()

setable_properties = ["box_l", "min_global_cut", "periodicity", "time",
                      "time_step", "timings"]
@highlander
class System(object):
    # NOTE: every attribute has to be declared at the class level.
    # This means that methods cannot define an attribute by using
    # `self.new_attr = somevalue` without declaring it inside this
    # indentation level, either as method, property or reference.
    part = particle_data.ParticleList()
    non_bonded_inter = interactions.NonBondedInteractions()
    bonded_inter = interactions.BondedInteractions()
    cell_system = CellSystem()
    thermostat = Thermostat()
    minimize_energy = MinimizeEnergy()
    actors = None
    analysis = None
    galilei = GalileiTransform()
    integrator = integrate.Integrator()
    if "CONSTRAINTS" in code_info.features():
        constraints = Constraints()
    lbboundaries = LBBoundaries()
    ekboundaries = EKBoundaries()

    auto_update_observables = AutoUpdateObservables()
    auto_update_correlators = AutoUpdateCorrelators()

    def __init__(self):
        self.actors = Actors(_system=self)
        self.analysis = Analysis(self)

    def __setattr__(self, name, value):
        if hasattr(self, name):
            object.__setattr__(self, name, value)
        else:
            raise AttributeError(
                "System does not have the attribute " + name + "."
                + "\nIf you know what you're doing, use "
                + "system.create_attr()"
            )

    def create_attr(self, *args, **kwargs):
        """
        Circumvents the __setattr__ lock. Allows to create and set new
        attributes to System instances.
        For *args, it initializes an attribute with None value, if the
        attribute does not already exist.
        For **kwargs, it simply calls super().__setattr__.
        """
        for arg in args:
            try: name = str(arg)
            except ValueError:
                print("Please pass either **kwargs or string *args to"
                      + "create_attr()"
                      )
                continue

            if hasattr(self, name):
                print("Attribute " + name + " already exists.")
            else:
                object.__setattr__(self, name, None)

        for name, value in list(kwargs.items()):
            object.__setattr__(self, name, value)

    # __getstate__ and __setstate__ define the pickle interaction
    def __getstate__(self):
        odict = {}
        for property_ in setable_properties:
            odict[property_] = System.__getattribute__(self, property_)
        return odict

    def __setstate__(self, params):
        for property_ in params.keys():
            System.__setattr__(self, property_, params[property_])

    @property
    def box_l(self):
        return np.array([box_l[0], box_l[1], box_l[2]])

    @box_l.setter
    def __set__(self, _box_l):
        if len(_box_l) != 3:
            raise ValueError("Box length must be of length 3")
        for i in range(3):
            if _box_l[i] <= 0:
                raise ValueError(
                    "Box length must be > 0 in all directions")
            box_l[i] = _box_l[i]

        mpi_bcast_parameter(0)

    @property
    def integ_switch(self):
        return integ_switch

    @property
    def periodicity(self, _periodic):
        global periodic
        if len(_periodic) != 3:
            raise ValueError(
                "periodicity must be of length 3, got length " + str(len(_periodic)))
        periodicity = np.zeros(3)
        for i in range(3):
            if _periodic[i] != 1:
                if "PARTIAL_PERIODIC" in code_info.features():
                    pass
                else:
                    raise ValueError(
                        "The feature PARTIAL_PERIODIC needs to be activated in myconfig.hpp")
        for i in range(3):
            periodicity[i] = _periodic[i]
        periodic = 4 * _periodic[2] + 2 * _periodic[1] + _periodic[0]
        # first 3 bits of periodic determine the periodicity
        mpi_bcast_parameter(FIELD_PERIODIC)

    @periodicity.setter
    def __get__(self):
        periodicity = np.zeros(3)
        periodicity[0] = periodic % 2
        periodicity[1] = int(periodic / 2) % 2
        periodicity[2] = int(periodic / 4) % 2
        return periodicity

    @property
    def time(self):
        global sim_time
        return sim_time

    @time.setter
    def __set__(self, _time):
        if _time < 0:
            raise ValueError("Simulation time must be >= 0")
        global sim_time
        sim_time = _time
        mpi_bcast_parameter(FIELD_SIMTIME)

    @property
    def smaller_time_step(self):
        return smaller_time_step

    @smaller_time_step.setter
    def __set__(self, _smaller_time_step):
        if "MULTI_TIMESTEP" in code_info.features():
            global smaller_time_step
            if _smaller_time_step <= 0:
                raise ValueError("Smaller time step must be positive")
            mpi_set_smaller_time_step(_smaller_time_step)

    @property
    def time_step(self):
        return time_step

    @time_step.setter
    def __set__(self, _time_step):
        if "LB" in code_info.features():
            global lbpar
        if "LB_GPU" in code_info.features():
            global lbpar_gpu
        if _time_step <= 0:
            raise ValueError("Time Step must be positive")
        if "LB" in code_info.features():
            if lbpar.tau >= 0.0 and _time_step < lbpar.tau:
                raise ValueError(
                    "Time Step (" + str(time_step) + ") must be > LB_time_step (" + str(lbpar.tau) + ")")
        if "LB_GPU" in code_info.features():
            if lbpar_gpu.tau >= 0.0 and _time_step < lbpar_gpu.tau:
                raise ValueError(
                    "Time Step (" + str(time_step) + ") must be > LB_time_step (" + str(lbpar_gpu.tau) + ")")
        mpi_set_time_step(_time_step)

    @property
    def timings(self):
        return timing_samples

    @timings.setter
    def __set__(self, _timings):
        global timing_samples
        if _timings <= 0:
            timing_samples = 0
        else:
            timing_samples = _timings

    @property
    def transfer_rate(self):
        return transfer_rate

    @property
    def max_cut_nonbonded(self):
        return max_cut_nonbonded

    @property
    def lattice_switch(self):
        return lattice_switch

    @property
    def max_cut_bonded(self):
        return max_cut_bonded

    @property
    def min_global_cut(self):
        return min_global_cut

    @min_global_cut.setter
    def __set__(self, _min_global_cut):
        global min_global_cut
        min_global_cut = _min_global_cut
        mpi_bcast_parameter(FIELD_MIN_GLOBAL_CUT)


    __seed = None

    def _get_PRNG_state_size(self):
        return get_state_size_of_generator()

    def set_random_state_PRNG(self):
        _state_size_plus_one = self._get_PRNG_state_size() + 1
        states = string_vec(n_nodes)
        rng = random.SystemRandom()  # true RNG that uses os.urandom()
        for i in range(n_nodes):
            states_on_node_i = []
            for j in range(_state_size_plus_one + 1):
                states_on_node_i.append(rng.randint(0, sys.maxint))
            states[i] = " ".join(map(str, states_on_node_i))
        mpi_random_set_stat(states)

    @property
    def seed(self):
        return self.__seed

    @seed.setter
    def __set__(self, _seed):
        seed_array = np.array([], dtype=int)
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

    @property
    def random_number_generator_state(self):
        rng_state = map(int, (mpi_random_get_stat().c_str()).split())
        return rng_state

    @random_number_generator_state.setter
    # sets the random number generator state in the core. this is of
    # interest for deterministic checkpointing
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


    def change_volume_and_rescale_particles(d_new, dir="xyz"):
        """Change box size and rescale particle coordinates
           change_volume_and_rescale_particles(d_new, dir="xyz")
           d_new: new length, dir=coordinate tow work on, "xyz" for isotropic
        """
        if d_new < 0:
            raise ValueError("No negative lengths")
        if dir == "xyz":
            d_new = d_new**(1. / 3.)
            rescale_boxl(3, d_new)
        elif dir == "x":
            rescale_boxl(0, d_new)
        elif dir == "y":
            rescale_boxl(1, d_new)
        elif dir == "z":
            rescale_boxl(2, d_new)
        else:
            raise ValueError(
                'Usage: changeVolume { <V_new> | <L_new> { "x" | "y" | "z" | "xyz" } }')

    def volume(self):
        """Return box volume"""
        return self.box_l[0] * self.box_l[1] * self.box_l[2]

    def distance(self, p1, p2):
        """Return the distance between the particles, respecting periodic boundaries"""
        res, a, b = np.zeros(3), np.zeros(3), np.zeros(3)
        a = p1.pos
        b = p2.pos
        get_mi_vector(res, a, b)
        return np.sqrt(res[0]**2 + res[1]**2 + res[2]**2)


if "CUDA" in code_info.features():
    cu = cuda_init.CudaInitHandle()
