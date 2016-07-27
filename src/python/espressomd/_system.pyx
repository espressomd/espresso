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
include "myconfig.pxi"

from globals cimport *
import numpy as np

cimport integrate
import integrate
import interactions
from actors import Actors
cimport cuda_init
import particle_data
import cuda_init
import code_info
from thermostat import Thermostat
from cellsystem import CellSystem
from minimize_energy import MinimizeEnergy
from polymer import Polymer
from analyze import Analysis
from galilei import GalileiTransform

import sys

setable_properties = ["box_l", "periodicity", "time",
                      "time_step", "timings"]

cdef class System:
    # NOTE: every attribute has to be declared at the class level.
    # This means that methods cannot define an attribute by using
    # `self.new_attr = somevalue` without declaring it inside this
    # indentation level, either as method, property or reference.
    doge = 1
    part = particle_data.ParticleList()
    non_bonded_inter = interactions.NonBondedInteractions()
    bonded_inter = interactions.BondedInteractions()
    cell_system = CellSystem()
    thermostat = Thermostat()
    minimize_energy = MinimizeEnergy()
    polymer = Polymer()
    actors = None
    analysis =None
    galilei = GalileiTransform()
    integrator = integrate.Integrator()

    def __init__(self):
        self.actors = Actors(_system=self)
        self.analysis=Analysis(self)

#        self.part = particle_data.particleList()
#        self.non_bonded_inter = interactions.NonBondedInteractions()
#        self.bonded_inter = interactions.BondedInteractions()

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
        def __set__(self, _box_l):
            if len(_box_l) != 3:
                raise ValueError("Box length must be of length 3")
            for i in range(3):
                if _box_l[i] <= 0:
                    raise ValueError(
                        "Box length must be > 0 in all directions")
                box_l[i] = _box_l[i]

            mpi_bcast_parameter(0)

        def __get__(self):
            return np.array([box_l[0], box_l[1], box_l[2]])

    property integ_switch:
        def __get__(self):
            return integ_switch

    property periodicity:
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
                            "No PARTIAL_PERIODIC activated")
            for i in range(3):
                periodicity[i] = _periodic[i]
            periodic = 4 * _periodic[2] + 2 * _periodic[1] + _periodic[0]
            # first 3 bits of periodic determine the periodicity
            # until we can handle contitional compilatio, periodic=7 is the
            # only value which makes sense
            mpi_bcast_parameter(FIELD_PERIODIC)

        def __get__(self):
            periodicity = np.zeros(3)
            periodicity[0] = periodic % 2
            periodicity[1] = int(periodic / 2) % 2
            periodicity[2] = int(periodic / 4) % 2
            return periodicity


    property time:
        def __set__(self, double _time):
            if _time <= 0:
                raise ValueError("Simulation time must be >= 0")
            global sim_time
            sim_time = _time
            mpi_bcast_parameter(FIELD_SIMTIME)

        def __get__(self):
            global sim_time
            return sim_time

    property smaller_time_step:
        def __set__(self, double _smaller_time_step):
            IF MULTI_TIMESTEP:
                global smaller_time_step
                if _smaller_time_step <= 0:
                    raise ValueError("Smaller time step must be positive")
                mpi_set_smaller_time_step(_smaller_time_step)

        def __get__(self):
            return smaller_time_step

    property time_step:
        def __set__(self, double _time_step):
            IF LB:
                global lbpar
            IF LB_GPU:
                global lbpar_gpu
            if _time_step <= 0:
                raise ValueError("Time Step must be positive")
            IF LB:
                if lbpar.tau >= 0.0 and _time_step > lbpar.tau:
                    raise ValueError(
                        "Time Step must be > LB_time_step (" + str(lbpar.tau) + ")")
            IF LB_GPU:
                if lbpar_gpu.tau >= 0.0 and _time_step > lbpar_gpu.tau:
                    raise ValueError(
                        "Time Step must be > LB_time_step (" + str(lbpar_gpu.tau) + ")")
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

    property transfer_rate:
        def __get__(self):
            return transfer_rate

    property max_cut_nonbonded:
        def __get__(self):
            return max_cut_nonbonded

    property lattice_switch:
        def __get__(self):
            return lattice_switch


    property max_cut_bonded:
        def __get__(self):
            return max_cut_bonded

    __seed = None
    property seed:
        def __set__(self, _seed):
            cdef vector[int] seed_array
            self.__seed = _seed
            if(isinstance(_seed, int) and self.n_nodes == 1):
                seed_array.resize(1)
                seed_array[0] = int(_seed)
                mpi_random_seed(0, seed_array)
            elif(hasattr(_seed, "__iter__")):
                if(len(_seed) < self.n_nodes or len(_seed) > self.n_nodes):
                    raise ValueError(
                        "The list needs to contain one seed value per node")
                seed_array.resize(len(_seed))
                for i in range(len(_seed)):
                    seed_array[i] = int(_seed[i])

                mpi_random_seed(self.n_nodes, seed_array)
            else:
                raise ValueError(
                    "The seed has to be an integer or a list of integers with one integer per node")

        def __get__(self):
            return self.__seed

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


# lbfluid=lb.DeviceList()
IF CUDA == 1:
    cu = cuda_init.CudaInitHandle()
