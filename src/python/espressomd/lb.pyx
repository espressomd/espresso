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
include "myconfig.pxi"
import os
import cython
import itertools
import numpy as np
cimport numpy as np
from libc cimport stdint
from .actors cimport Actor
from . cimport cuda_init
from . import cuda_init
from . import utils
from .utils import array_locked, is_valid_type
from .utils cimport Vector3i, Vector3d, Vector6d, Vector19d, make_array_locked
from .globals cimport time_step


def _construct(cls, params):
    obj = cls(**params)
    obj._params = params
    return obj


def assert_agrid_tau_set(obj):
    assert obj.agrid != obj.default_params()['agrid'] and obj.tau != obj.default_params()[
        'tau'], "tau and agrid have to be set first!"


cdef class HydrodynamicInteraction(Actor):
    def _lb_init(self):
        raise Exception(
            "Subclasses of HydrodynamicInteraction must define the _lb_init() method.")

    def __reduce__(self):
        return _construct, (self.__class__, self._params), None

    def __getitem__(self, key):
        if isinstance(key, (tuple, list, np.ndarray)):
            if len(key) == 3:
                return LBFluidRoutines(np.array(key))
        else:
            raise Exception(
                "%s is not a valid key. Should be a point on the nodegrid e.g. lbf[0,0,0]," % key)

    # validate the given parameters on actor initialization
    ####################################################
    def validate_params(self):
        default_params = self.default_params()

        utils.check_type_or_throw_except(
            self._params["kT"], 1, float, "kT must be a number")
        if self._params["kT"] > 0. and not self._params["seed"]:
            raise ValueError(
                "seed has to be given if temperature is not 0.")

        if self._params["dens"] == default_params["dens"]:
            raise Exception("LB_FLUID density not set")
        elif not (self._params["dens"] > 0.0 and (is_valid_type(self._params["dens"], float) or is_valid_type(self._params["dens"], int))):
            raise ValueError("Density must be a positive double")

        if self._params["tau"] <= 0.:
            raise ValueError("tau has to be a positive double")

    def valid_keys(self):
        return "agrid", "dens", "ext_force_density", "visc", "tau", "bulk_visc", "gamma_odd", "gamma_even", "kT", "seed"

    def required_keys(self):
        return ["dens", "agrid", "visc", "tau"]

    def default_params(self):
        return {"agrid": -1.0,
                "dens": -1.0,
                "ext_force_density": [0.0, 0.0, 0.0],
                "visc": -1.0,
                "bulk_visc": -1.0,
                "tau": -1.0,
                "seed": None,
                "kT": 0.}

    def _set_lattice_switch(self):
        raise Exception(
            "Subclasses of HydrodynamicInteraction must define the _set_lattice_switch() method.")

    def _set_params_in_es_core(self):
        default_params = self.default_params()
        self.agrid = self._params['agrid']
        self.tau = self._params['tau']
        self.density = self._params['dens']

        if self._params['kT'] > 0.:
            self.seed = self._params['seed']
        self.kT = self._params['kT']

        self.viscosity = self._params['visc']
        if self._params['bulk_visc'] != default_params['bulk_visc']:
            self.bulk_viscosity = self._params['bulk_visc']

        self.ext_force_density = self._params["ext_force_density"]

        if "gamma_odd" in self._params:
            python_lbfluid_set_gamma_odd(self._params["gamma_odd"])

        if "gamma_even" in self._params:
            python_lbfluid_set_gamma_even(self._params["gamma_even"])

        lb_lbfluid_sanity_checks()
        utils.handle_errors("LB fluid activation")

    def _get_params_from_es_core(self):
        default_params = self.default_params()
        self._params['agrid'] = self.agrid
        self._params["tau"] = self.tau
        self._params['dens'] = self.density
        self._params["kT"] = self.kT
        if self._params['kT'] > 0.0:
            self._params['seed'] = self.seed
        self._params['visc'] = self.viscosity
        if not self._params["bulk_visc"] == default_params["bulk_visc"]:
            self._params['bulk_visc'] = self.bulk_viscosity
        self._params['ext_force_density'] = self.ext_force_density

        return self._params

    def set_interpolation_order(self, interpolation_order):
        """ Set the order for the fluid interpolation scheme.

        Parameters
        ----------
        interpolation_order : :obj:`str`
            ``"linear"`` for linear interpolation, ``"quadratic"`` for quadratic interpolation.

        """
        if (interpolation_order == "linear"):
            lb_lbinterpolation_set_interpolation_order(linear)
        elif (interpolation_order == "quadratic"):
            lb_lbinterpolation_set_interpolation_order(quadratic)
        else:
            raise ValueError("Invalid parameter")

    def get_interpolated_velocity(self, pos):
        """Get LB fluid velocity at specified position.

        Parameters
        ----------
        pos : (3,) array_like of :obj:`float`
              The position at which velocity is requested.

        Returns
        -------
        v : (3,) array_like :obj:`float`
            The LB fluid velocity at ``pos``.

        """
        cdef Vector3d p

        for i in range(3):
            p[i] = pos[i]
        cdef Vector3d v = lb_lbfluid_get_interpolated_velocity(p) * lb_lbfluid_get_lattice_speed()
        return make_array_locked(v)

    def print_vtk_velocity(self, path, bb1=None, bb2=None):
        cdef vector[int] bb1_vec
        cdef vector[int] bb2_vec
        if bb1 is None or bb2 is None:
            lb_lbfluid_print_vtk_velocity(utils.to_char_pointer(path))
        else:
            bb1_vec = bb1
            bb2_vec = bb2
            lb_lbfluid_print_vtk_velocity(
                utils.to_char_pointer(path), bb1_vec, bb2_vec)

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

    def _activate_method(self):
        raise Exception(
            "Subclasses of HydrodynamicInteraction have to implement _activate_method.")

    def _deactivate_method(self):
        lb_lbfluid_set_lattice_switch(NONE)

    property shape:
        def __get__(self):
            cdef Vector3i shape = lb_lbfluid_get_shape()
            return (shape[0], shape[1], shape[2])

    property kT:
        def __get__(self):
            return lb_lbfluid_get_kT()

        def __set__(self, kT):
            cdef double _kT = kT
            lb_lbfluid_set_kT(_kT)

    property seed:
        def __get__(self):
            return lb_lbfluid_get_rng_state()

        def __set__(self, seed):
            cdef stdint.uint64_t _seed = seed
            lb_lbfluid_set_rng_state(seed)

    property stress:
        def __get__(self):
            cdef Vector6d stress = python_lbfluid_get_stress(self.agrid, self.tau)
            return array_locked(np.array([[stress[0], stress[1], stress[3]],
                                          [stress[1], stress[2], stress[4]],
                                          [stress[3], stress[4], stress[5]]]))

        def __set__(self, value):
            raise NotImplementedError

    property ext_force_density:
        def __get__(self):
            assert_agrid_tau_set(self)
            cdef Vector3d res
            res = python_lbfluid_get_ext_force_density(
                self.agrid, self.tau)
            return make_array_locked(res)

        def __set__(self, ext_force_density):
            assert_agrid_tau_set(self)
            python_lbfluid_set_ext_force_density(
                ext_force_density, self.agrid, self.tau)

    property density:
        def __get__(self):
            assert_agrid_tau_set(self)
            return python_lbfluid_get_density(self.agrid)

        def __set__(self, density):
            assert_agrid_tau_set(self)
            python_lbfluid_set_density(density, self.agrid)

    property viscosity:
        def __get__(self):
            assert_agrid_tau_set(self)
            return python_lbfluid_get_viscosity(self.agrid, self.tau)

        def __set__(self, viscosity):
            assert_agrid_tau_set(self)
            python_lbfluid_set_viscosity(viscosity, self.agrid, self.tau)

    property bulk_viscosity:
        def __get__(self):
            assert_agrid_tau_set(self)
            return python_lbfluid_get_bulk_viscosity(self.agrid, self.tau)

        def __set__(self, viscosity):
            assert_agrid_tau_set(self)
            python_lbfluid_set_bulk_viscosity(viscosity, self.agrid, self.tau)

    property tau:
        def __get__(self):
            return lb_lbfluid_get_tau()

        def __set__(self, tau):
            lb_lbfluid_set_tau(tau)
            if time_step > 0.0:
                check_tau_time_step_consistency(tau, time_step)

    property agrid:
        def __get__(self):
            return lb_lbfluid_get_agrid()

        def __set__(self, agrid):
            lb_lbfluid_set_agrid(agrid)

    def nodes(self):
        """Provides a generator for iterating over all lb nodes"""

        shape = self.shape
        for i, j, k in itertools.product(
                range(shape[0]), range(shape[1]), range(shape[2])):
            yield self[i, j, k]


cdef class LBFluid(HydrodynamicInteraction):
    """
    Initialize the lattice-Boltzmann method for hydrodynamic flow using the CPU.

    """

    def _set_lattice_switch(self):
        lb_lbfluid_set_lattice_switch(CPU)

    def _activate_method(self):
        self.validate_params()
        self._set_lattice_switch()
        self._set_params_in_es_core()

IF CUDA:
    cdef class LBFluidGPU(HydrodynamicInteraction):
        """
        Initialize the lattice-Boltzmann method for hydrodynamic flow using the GPU.

        """

        def _set_lattice_switch(self):
            lb_lbfluid_set_lattice_switch(GPU)

        def _activate_method(self):
            self.validate_params()
            self._set_lattice_switch()
            self._set_params_in_es_core()

        @cython.boundscheck(False)
        @cython.wraparound(False)
        def get_interpolated_fluid_velocity_at_positions(self, np.ndarray[double, ndim=2, mode="c"] positions not None, three_point=False):
            """Calculate the fluid velocity at given positions.

            Parameters
            ----------
            positions : (N,3) numpy-array of type :obj:`float`
                The 3-dimensional positions.

            Returns
            -------
            velocities : (N,3) numpy-array of type :obj:`float`
                The 3-dimensional LB fluid velocities.

            Raises
            ------
            AssertionError
                If shape of ``positions`` not (N,3).

            """
            assert positions.shape[1] == 3, \
                "The input array must have shape (N,3)"
            cdef int length
            length = positions.shape[0]
            velocities = np.empty_like(positions)
            if three_point:
                quadratic_velocity_interpolation( < double * >np.PyArray_GETPTR2(positions, 0, 0), < double * >np.PyArray_GETPTR2(velocities, 0, 0), length)
            else:
                linear_velocity_interpolation( < double * >np.PyArray_GETPTR2(positions, 0, 0), < double * >np.PyArray_GETPTR2(velocities, 0, 0), length)
            return velocities * lb_lbfluid_get_lattice_speed()

cdef class LBFluidRoutines:
    cdef Vector3i node

    def __init__(self, key):
        utils.check_type_or_throw_except(
            key, 3, int, "The index of an lb fluid node consists of three integers.")
        self.node[0] = key[0]
        self.node[1] = key[1]
        self.node[2] = key[2]
        if not lb_lbnode_is_index_valid(self.node):
            raise ValueError("LB node index out of bounds")

    property index:
        def __get__(self):
            return (self.node[0], self.node[1], self.node[2])

    property velocity:
        def __get__(self):
            return make_array_locked(python_lbnode_get_velocity(self.node))

        def __set__(self, value):
            cdef Vector3d c_velocity
            if all(is_valid_type(v, float) for v in value) and len(value) == 3:
                c_velocity[0] = value[0]
                c_velocity[1] = value[1]
                c_velocity[2] = value[2]
                python_lbnode_set_velocity(self.node, c_velocity)
            else:
                raise ValueError(
                    "Velocity has to be of shape 3 and type float.")
    property density:
        def __get__(self):
            return python_lbnode_get_density(self.node)

        def __set__(self, value):
            python_lbnode_set_density(self.node, value)

    property stress:
        def __get__(self):
            cdef Vector6d stress = python_lbnode_get_stress(self.node)
            return array_locked(np.array([[stress[0], stress[1], stress[3]],
                                          [stress[1], stress[2], stress[4]],
                                          [stress[3], stress[4], stress[5]]]))

        def __set__(self, value):
            raise NotImplementedError

    property stress_neq:
        def __get__(self):
            cdef Vector6d stress = python_lbnode_get_stress_neq(self.node)
            return array_locked(np.array([[stress[0], stress[1], stress[3]],
                                          [stress[1], stress[2], stress[4]],
                                          [stress[3], stress[4], stress[5]]]))

        def __set__(self, value):
            raise NotImplementedError

    property population:
        def __get__(self):
            cdef Vector19d double_return
            double_return = lb_lbnode_get_pop(self.node)
            return array_locked(np.array([double_return[0],
                                          double_return[1],
                                          double_return[2],
                                          double_return[3],
                                          double_return[4],
                                          double_return[5],
                                          double_return[6],
                                          double_return[7],
                                          double_return[8],
                                          double_return[9],
                                          double_return[10],
                                          double_return[11],
                                          double_return[12],
                                          double_return[13],
                                          double_return[14],
                                          double_return[15],
                                          double_return[16],
                                          double_return[17],
                                          double_return[18]]
                                         ))

        def __set__(self, population):
            cdef Vector19d _population
            for i in range(19):
                _population[i] = population[i]
            lb_lbnode_set_pop(self.node, _population)

    property boundary:
        def __get__(self):
            return lb_lbnode_get_boundary(self.node)

        def __set__(self, value):
            raise NotImplementedError
