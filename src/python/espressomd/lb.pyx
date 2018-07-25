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
from __future__ import print_function, absolute_import, division
include "myconfig.pxi"
import os
import cython
import numpy as np
cimport numpy as np
from .actors cimport Actor
from . cimport cuda_init
from . import cuda_init
from globals cimport *
from copy import deepcopy
from . import utils
from espressomd.utils import array_locked, is_valid_type

# Actor class
####################################################
cdef class HydrodynamicInteraction(Actor):
    def _lb_init(self):
        raise Exception(
            "Subclasses of HydrodynamicInteraction must define the _lb_init() method.")

# LBFluid main class
####################################################
IF LB_GPU or LB:
    cdef class LBFluid(HydrodynamicInteraction):
        """
        Initialize the lattice-Boltzmann method for hydrodynamic flow using the CPU.

        """            

        def __getitem__(self, key):
            if isinstance(key, tuple) or isinstance(key, list) or isinstance(key, np.ndarray):
                if len(key) == 3:
                    return LBFluidRoutines(np.array(key))
            else: 
                raise Exception("%s is not a valid key. Should be a point on the nodegrid e.g. lbf[0,0,0]," %key)


        # validate the given parameters on actor initalization
        ####################################################
        def validate_params(self):
            default_params = self.default_params()

            IF SHANCHEN:
                if not hasattr(self._params["dens"],"__getitem__"):
                    raise ValueError(
                        "Density must be two positive double (ShanChen)")
                if not (self._params["dens"][0] > 0.0 and self._params["dens"][1] > 0.0):
                    raise ValueError(
                        "Density must be two positive double (ShanChen)")
            ELSE:
                if self._params["dens"] == default_params["dens"]:
                    raise Exception("LB_FLUID density not set")
                else:
                    if not (self._params["dens"] > 0.0 and (is_valid_type(self._params["dens"], float) or is_valid_type(self._params["dens"], int))):
                        raise ValueError("Density must be one positive double")

        # list of valid keys for parameters
        ####################################################
        def valid_keys(self):
            return "agrid", "dens", "fric", "ext_force_density", "visc", "tau", "couple"

        # list of esential keys required for the fluid
        ####################################################
        def required_keys(self):
            return ["dens", "agrid", "visc", "tau"]

        # list of default parameters
        ####################################################
        def default_params(self):
            IF SHANCHEN:
                return {"agrid": -1.0,
                        "dens": [-1.0, -1.0],
                        "fric": [-1.0, -1.0],
                        "ext_force_density": [0.0, 0.0, 0.0],
                        "visc": [-1.0, -1.0],
                        "bulk_visc": [-1.0, -1.0],
                        "tau": -1.0,
                        "couple": "2pt"}
            ELSE:
                return {"agrid": -1.0,
                        "dens": -1.0,
                        "fric": -1.0,
                        "ext_force_density": [0.0, 0.0, 0.0],
                        "visc": -1.0,
                        "bulk_visc": -1.0,
                        "tau": -1.0,
                        "couple": "2pt"}

        # function that calls wrapper functions which set the parameters at C-Level
        ####################################################
        def _set_lattice_switch(self):
            if lb_set_lattice_switch(1):
                raise Exception("lb_set_lattice_switch error")

        def _set_params_in_es_core(self):
            default_params = self.default_params()

            if python_lbfluid_set_density(self._params["dens"]):
                raise Exception("lb_lbfluid_set_density error")

            if python_lbfluid_set_tau(self._params["tau"]):
                raise Exception("lb_lbfluid_set_tau error")

            if python_lbfluid_set_visc(self._params["visc"]):
                raise Exception("lb_lbfluid_set_visc error")

            if self._params["bulk_visc"] != self.default_params()["bulk_visc"]:
                if python_lbfluid_set_bulk_visc(self._params["bulk_visc"]):
                    raise Exception("lb_lbfluid_set_bulk_visc error")

            if python_lbfluid_set_agrid(self._params["agrid"]):
                raise Exception("lb_lbfluid_set_agrid error")

            if self._params["fric"]!=default_params["fric"]:
                if python_lbfluid_set_friction(self._params["fric"]):
                    raise Exception("lb_lbfluid_set_friction error")

            if python_lbfluid_set_ext_force_density(self._params["ext_force_density"]):
                raise Exception("lb_lbfluid_set_ext_force_density error")

            if python_lbfluid_set_couple_flag(self._params["couple"]):
                raise Exception("lb_lbfluid_set_couple_flag error")
            utils.handle_errors("LB fluid activation")

        # function that calls wrapper functions which get the parameters from C-Level
        ####################################################
        def _get_params_from_es_core(self):
            default_params = self.default_params()

            if python_lbfluid_get_density(self._params["dens"]):
                raise Exception("lb_lbfluid_get_density error")

            if python_lbfluid_get_tau(self._params["tau"]):
                raise Exception("lb_lbfluid_set_tau error")

            if python_lbfluid_get_visc(self._params["visc"]):
                raise Exception("lb_lbfluid_set_visc error")

            if not self._params["bulk_visc"] == default_params["bulk_visc"]:
                if python_lbfluid_get_bulk_visc(self._params["bulk_visc"]):
                    raise Exception("lb_lbfluid_set_bulk_visc error")

            if python_lbfluid_get_agrid(self._params["agrid"]):
                raise Exception("lb_lbfluid_set_agrid error")

            if not self._params["fric"] == default_params["fric"]:
                if python_lbfluid_get_friction(self._params["fric"]):
                    raise Exception("lb_lbfluid_set_friction error")

            if not self._params["ext_force_density"] == default_params["ext_force_density"]:
                if python_lbfluid_get_ext_force_density(self._params["ext_force_density"]):
                    raise Exception("lb_lbfluid_set_ext_force_density error")

            if not self._params["couple"] == default_params["couple"]:
                if python_lbfluid_get_couple_flag(self._params["couple"]):
                    raise Exception("lb_lbfluid_get_couple_flag error")

            return self._params
        
        def get_interpolated_velocity(self, pos):
            """Get LB fluid velocity at specified position.

            Parameters
            ----------
            pos : array_like :obj:`float`
                  The position at which velocity is requested.

            Returns
            -------
            v : array_like :obj:`float`
                The LB fluid velocity at ``pos``.

            """
            cdef Vector3d p	    
            cdef double[3] v

            for i in range(3):
                p[i] = pos[i]

            lb_lbfluid_get_interpolated_velocity_global(p, v)
            return v

        # input/output function wrappers for whole LB fields
        ####################################################
        def print_vtk_velocity(self, path, bb1=None, bb2=None):
            cdef vector[int] bb1_vec
            cdef vector[int] bb2_vec
            if bb1 is None or bb2 is None:
                lb_lbfluid_print_vtk_velocity(utils.to_char_pointer(path))
            else:
                bb1_vec = bb1
                bb2_vec = bb2
                lb_lbfluid_print_vtk_velocity(utils.to_char_pointer(path), bb1_vec, bb2_vec)
        def print_vtk_boundary(self, path):
            lb_lbfluid_print_vtk_boundary(utils.to_char_pointer(path))
        def print_velocity(self, path):
            lb_lbfluid_print_velocity(utils.to_char_pointer(path))
        def print_boundary(self, path):
            lb_lbfluid_print_boundary(utils.to_char_pointer(path))
        def save_checkpoint(self, path, binary):
            tmp_path = path + ".__tmp__"
            lb_lbfluid_save_checkpoint(utils.to_char_pointer(tmp_path), binary)
            os.rename(tmp_path, path)
        def load_checkpoint(self, path, binary):
            lb_lbfluid_load_checkpoint(utils.to_char_pointer(path), binary)

        # Activate Actor
        ####################################################
        def _activate_method(self):
            self.validate_params()
            self._set_lattice_switch()
            self._set_params_in_es_core()
            utils.handle_errors("LB fluid activation")
            IF LB:
                return
            ELSE:
                raise Exception("LB not compiled in")

        def _deactivate_method(self):
            if lb_set_lattice_switch(0):
                raise Exception("lb_set_lattice_switch error")

IF LB_GPU:
    cdef class LBFluidGPU(LBFluid):
        """
        Initialize the lattice-Boltzmann method for hydrodynamic flow using the GPU.

        """

        def remove_total_momentum(self):
            lb_lbfluid_remove_total_momentum()
        
        def _set_lattice_switch(self):
            if lb_set_lattice_switch(2):
                raise Exception("lb_set_lattice_switch error")


        def _activate_method(self):
            self.validate_params()
            self._set_lattice_switch()
            self._set_params_in_es_core()
            IF LB_GPU:
                return
            ELSE:
                raise Exception("LB_GPU not compiled in")


        @cython.boundscheck(False)
        @cython.wraparound(False)
        def get_interpolated_fluid_velocity_at_positions(self, np.ndarray[double, ndim=2, mode="c"] positions not None):
            """Calculate the fluid velocity at given positions.

            Parameters
            ----------
            positions : numpy-array of type :obj:`float` of shape (N,3)
                        The 3-dimensional positions.

            Returns
            -------
            velocities : numpy-array of type :obj:`float` of shape (N,3)
                         The 3-dimensional LB fluid velocities.

            Raises
            ------
            AssertionError
                If shape of ``positions`` not (N,3).

            """
            assert positions.shape[1] == 3, "The input array must have shape (N,3)"
            cdef int length
            length = positions.shape[0]
            velocities = np.empty_like(positions)
            lb_lbfluid_get_interpolated_velocity_at_positions(<double *>np.PyArray_GETPTR2(positions, 0, 0), <double *>np.PyArray_GETPTR2(velocities, 0, 0), length)
            return velocities

IF LB or LB_GPU:
    cdef class LBFluidRoutines(object):
        cdef int node[3]
        def __init__(self, key):
            utils.check_type_or_throw_except(key,3,int, "The index of an lb fluid node consists of three integers.")
            self.node[0] = key[0]
            self.node[1] = key[1]
            self.node[2] = key[2]

        property velocity:
            def __get__(self):
                cdef double[3] double_return
                lb_lbnode_get_u(self.node, double_return)
                return array_locked(double_return)

            def __set__(self, value):
                cdef double[3] host_velocity
                if all(is_valid_type(v, float) for v in value) and len(value) == 3:
                    host_velocity = value
                    lb_lbnode_set_u(self.node, host_velocity)
                else:
                    raise ValueError("Velocity has to be of shape 3 and type float.")
        property density:
            def __get__(self):
                cdef double[3] double_return
                lb_lbnode_get_rho(self.node, double_return)
                return array_locked(double_return)

            def __set__(self, value):
                raise NotImplementedError


        property pi:
            def __get__(self):
                cdef double[6] pi
                lb_lbnode_get_pi(self.node, pi)
                return array_locked(np.array([[pi[0],pi[1],pi[3]],
                                              [pi[1],pi[2],pi[4]],
                                              [pi[3],pi[4],pi[5]]]))

            def __set__(self, value):
                raise NotImplementedError

        property pi_neq:
            def __get__(self):
                cdef double[6] pi
                lb_lbnode_get_pi_neq(self.node, pi)
                return array_locked(np.array([[pi[0],pi[1],pi[3]],
                                              [pi[1],pi[2],pi[4]],
                                              [pi[3],pi[4],pi[5]]]))

            def __set__(self, value):
                raise NotImplementedError

        property population:
            def __get__(self):
                cdef double[19] double_return
                lb_lbnode_get_pop(self.node, double_return)
                return array_locked(double_return)

            def __set__(self, value):
                cdef double[19] double_return = value
                lb_lbnode_set_pop(self.node, double_return)

        property boundary:
            def __get__(self):
                cdef int int_return
                lb_lbnode_get_boundary(self.node, &int_return)
                return int_return

            def __set__(self, value):
                raise NotImplementedError
