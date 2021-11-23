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
import functools
import numpy as np
cimport numpy as np
from libc cimport stdint
from .actors cimport Actor
from . cimport cuda_init
from . import cuda_init
from . import utils
from .utils import array_locked, is_valid_type, check_type_or_throw_except
from .utils cimport Vector3i, Vector3d, Vector6d, Vector19d, make_array_locked, make_Vector3d
from .integrate cimport get_time_step


def _construct(cls, params):
    obj = cls(**params)
    obj._params = params
    return obj


cdef class HydrodynamicInteraction(Actor):
    """
    Base class for LB implementations.

    Parameters
    ----------
    agrid : :obj:`float`
        Lattice constant. The box size in every direction must be an integer
        multiple of ``agrid``.
    tau : :obj:`float`
        LB time step, must be an integer multiple of the MD time step.
    dens : :obj:`float`
        Fluid density.
    visc : :obj:`float`
        Fluid kinematic viscosity.
    bulk_visc : :obj:`float`, optional
        Fluid bulk viscosity.
    gamma_odd : :obj:`int`, optional
        Relaxation parameter :math:`\\gamma_{\\textrm{odd}}` for kinetic modes.
    gamma_even : :obj:`int`, optional
        Relaxation parameter :math:`\\gamma_{\\textrm{even}}` for kinetic modes.
    ext_force_density : (3,) array_like of :obj:`float`, optional
        Force density applied on the fluid.
    kT : :obj:`float`, optional
        Thermal energy of the simulated heat bath (for thermalized fluids).
        Set it to 0 for an unthermalized fluid.
    seed : :obj:`int`, optional
        Initial counter value (or seed) of the philox RNG.
        Required for a thermalized fluid. Must be positive.
    """

    def _assert_agrid_tau_set(self):
        unset = self.default_params()
        assert self.agrid != unset['agrid'] and self.tau != unset['tau'], \
            "tau and agrid have to be set first!"

    def _lb_init(self):
        raise Exception(
            "Subclasses of HydrodynamicInteraction must define the _lb_init() method.")

    def __reduce__(self):
        return _construct, (self.__class__, self._params), None

    def __getitem__(self, key):
        cdef Vector3i shape
        if isinstance(key, (tuple, list, np.ndarray)):
            if len(key) == 3:
                if any(isinstance(typ, slice) for typ in key):
                    shape = lb_lbfluid_get_shape()
                    return LBSlice(key, (shape[0], shape[1], shape[2]))
                else:
                    return LBFluidRoutines(np.array(key))
        else:
            raise Exception(
                "%s is not a valid key. Should be a point on the nodegrid e.g. lbf[0,0,0], or a slice" % key)
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
        return {"agrid", "dens", "ext_force_density", "visc", "tau",
                "bulk_visc", "gamma_odd", "gamma_even", "kT", "seed"}

    def required_keys(self):
        return {"dens", "agrid", "visc", "tau"}

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
        if 'gamma_odd' in self._params:
            self._params['gamma_odd'] = lb_lbfluid_get_gamma_odd()
        if 'gamma_even' in self._params:
            self._params['gamma_even'] = lb_lbfluid_get_gamma_even()

        return self._params

    def set_interpolation_order(self, interpolation_order):
        """ Set the order for the fluid interpolation scheme.

        Parameters
        ----------
        interpolation_order : :obj:`str`, \{"linear", "quadratic"\}
            ``"linear"`` for trilinear interpolation, ``"quadratic"`` for
            quadratic interpolation. For the CPU implementation of LB, only
            ``"linear"`` is available.

        """
        if interpolation_order == "linear":
            lb_lbinterpolation_set_interpolation_order(linear)
        elif interpolation_order == "quadratic":
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
        return python_lbnode_get_interpolated_velocity(make_Vector3d(pos))

    def write_vtk_velocity(self, path, bb1=None, bb2=None):
        """Write the LB fluid velocity to a VTK file.
        If both ``bb1`` and ``bb2`` are specified, return a subset of the grid.

        Parameters
        ----------
        path : :obj:`str`
            Path to the output ASCII file.
        bb1 : (3,) array_like of :obj:`int`, optional
            Node indices of the lower corner of the bounding box.
        bb2 : (3,) array_like of :obj:`int`, optional
            Node indices of the upper corner of the bounding box.

        """
        cdef vector[int] bb1_vec
        cdef vector[int] bb2_vec
        if bb1 is None and bb2 is None:
            lb_lbfluid_print_vtk_velocity(utils.to_char_pointer(path))
        elif bb1 is None or bb2 is None:
            raise ValueError(
                "Invalid parameter: must provide either both bb1 and bb2, or none of them")
        else:
            check_type_or_throw_except(bb1, 3, int,
                                       "bb1 has to be an integer list of length 3")
            check_type_or_throw_except(bb2, 3, int,
                                       "bb2 has to be an integer list of length 3")
            bb1_vec = bb1
            bb2_vec = bb2
            lb_lbfluid_print_vtk_velocity(
                utils.to_char_pointer(path), bb1_vec, bb2_vec)

    def write_vtk_boundary(self, path):
        """Write the LB boundaries to a VTK file.

        Parameters
        ----------
        path : :obj:`str`
            Path to the output ASCII file.

        """
        lb_lbfluid_print_vtk_boundary(utils.to_char_pointer(path))

    def write_velocity(self, path):
        """Write the LB fluid velocity to a data file that can be loaded by
        numpy, with format "x y z vx vy vz".

        Parameters
        ----------
        path : :obj:`str`
            Path to the output data file.

        """
        lb_lbfluid_print_velocity(utils.to_char_pointer(path))

    def write_boundary(self, path):
        """Write the LB boundaries to a data file that can be loaded by numpy,
        with format "x y z u".

        Parameters
        ----------
        path : :obj:`str`
            Path to the output data file.

        """
        lb_lbfluid_print_boundary(utils.to_char_pointer(path))

    def save_checkpoint(self, path, binary):
        '''
        Write LB node populations to a file.
        :class:`~espressomd.lbboundaries.LBBoundaries`
        information is not written to the file.
        '''
        tmp_path = path + ".__tmp__"
        lb_lbfluid_save_checkpoint(utils.to_char_pointer(tmp_path), binary)
        os.rename(tmp_path, path)

    def load_checkpoint(self, path, binary):
        '''
        Load LB node populations from a file.
        :class:`~espressomd.lbboundaries.LBBoundaries`
        information is not available in the file. The boundary
        information of the grid will be set to zero,
        even if :class:`~espressomd.lbboundaries.LBBoundaries`
        contains :class:`~espressomd.lbboundaries.LBBoundary`
        objects (they are ignored).
        '''
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

    property pressure_tensor:
        def __get__(self):
            tensor = python_lbfluid_get_pressure_tensor(self.agrid, self.tau)
            return array_locked(tensor)

        def __set__(self, value):
            raise NotImplementedError

    property ext_force_density:
        def __get__(self):
            self._assert_agrid_tau_set()
            return python_lbfluid_get_ext_force_density(self.agrid, self.tau)

        def __set__(self, ext_force_density):
            self._assert_agrid_tau_set()
            python_lbfluid_set_ext_force_density(
                make_Vector3d(ext_force_density), self.agrid, self.tau)

    property density:
        def __get__(self):
            self._assert_agrid_tau_set()
            return python_lbfluid_get_density(self.agrid)

        def __set__(self, density):
            self._assert_agrid_tau_set()
            python_lbfluid_set_density(density, self.agrid)

    property viscosity:
        def __get__(self):
            self._assert_agrid_tau_set()
            return python_lbfluid_get_viscosity(self.agrid, self.tau)

        def __set__(self, viscosity):
            self._assert_agrid_tau_set()
            python_lbfluid_set_viscosity(viscosity, self.agrid, self.tau)

    property bulk_viscosity:
        def __get__(self):
            self._assert_agrid_tau_set()
            return python_lbfluid_get_bulk_viscosity(self.agrid, self.tau)

        def __set__(self, viscosity):
            self._assert_agrid_tau_set()
            python_lbfluid_set_bulk_viscosity(viscosity, self.agrid, self.tau)

    property tau:
        def __get__(self):
            return lb_lbfluid_get_tau()

        def __set__(self, tau):
            lb_lbfluid_set_tau(tau)
            if get_time_step() > 0.0:
                check_tau_time_step_consistency(tau, get_time_step())

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
    See :class:`HydrodynamicInteraction` for the list of parameters.

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
        See :class:`HydrodynamicInteraction` for the list of parameters.

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
                quadratic_velocity_interpolation(< double * >np.PyArray_GETPTR2(positions, 0, 0), < double * >np.PyArray_GETPTR2(velocities, 0, 0), length)
            else:
                linear_velocity_interpolation(< double * >np.PyArray_GETPTR2(positions, 0, 0), < double * >np.PyArray_GETPTR2(velocities, 0, 0), length)
            return velocities * lb_lbfluid_get_lattice_speed()

ELSE:
    cdef class LBFluidGPU(HydrodynamicInteraction):
        def __init__(self, *args, **kwargs):
            raise Exception("LBFluidGPU not compiled in.")


cdef class LBFluidRoutines:

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
            return python_lbnode_get_velocity(self.node)

        def __set__(self, value):
            utils.check_type_or_throw_except(
                value, 3, float, "velocity has to be 3 floats")
            python_lbnode_set_velocity(self.node, make_Vector3d(value))

    property density:
        def __get__(self):
            return python_lbnode_get_density(self.node)

        def __set__(self, value):
            python_lbnode_set_density(self.node, value)

    property pressure_tensor:
        def __get__(self):
            tensor = python_lbnode_get_pressure_tensor(self.node)
            return array_locked(tensor)

        def __set__(self, value):
            raise NotImplementedError

    property pressure_tensor_neq:
        def __get__(self):
            tensor = python_lbnode_get_pressure_tensor_neq(self.node)
            return array_locked(tensor)

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

    def __eq__(self, obj1):
        index_1 = np.array(self.index)
        index_2 = np.array(obj1.index)
        return all(index_1 == index_2)

    def __hash__(self):
        return hash(self.index)


class LBSlice:

    def __init__(self, key, shape):
        self.x_indices, self.y_indices, self.z_indices = self.get_indices(
            key, shape[0], shape[1], shape[2])

    def get_indices(self, key, shape_x, shape_y, shape_z):
        x_indices = np.atleast_1d(np.arange(shape_x)[key[0]])
        y_indices = np.atleast_1d(np.arange(shape_y)[key[1]])
        z_indices = np.atleast_1d(np.arange(shape_z)[key[2]])
        return x_indices, y_indices, z_indices

    def get_values(self, x_indices, y_indices, z_indices, prop_name):
        shape_res = np.shape(
            getattr(LBFluidRoutines(np.array([0, 0, 0])), prop_name))
        res = np.zeros(
            (x_indices.size,
             y_indices.size,
             z_indices.size,
             *shape_res))
        for i, x in enumerate(x_indices):
            for j, y in enumerate(y_indices):
                for k, z in enumerate(z_indices):
                    res[i, j, k] = getattr(LBFluidRoutines(
                        np.array([x, y, z])), prop_name)
        if shape_res == (1,):
            res = np.squeeze(res, axis=-1)
        return array_locked(res)

    def set_values(self, x_indices, y_indices, z_indices, prop_name, value):
        for i, x in enumerate(x_indices):
            for j, y in enumerate(y_indices):
                for k, z in enumerate(z_indices):
                    setattr(LBFluidRoutines(
                        np.array([x, y, z])), prop_name, value[i, j, k])

    def __iter__(self):
        indices = [(x, y, z) for (x, y, z) in itertools.product(
            self.x_indices, self.y_indices, self.z_indices)]
        return (LBFluidRoutines(np.array(index)) for index in indices)


def _add_lb_slice_properties():
    """
    Automatically add all of LBFluidRoutines's properties to LBSlice.

    """

    def set_attribute(lb_slice, value, attribute):
        """
        Setter function that sets attribute on every member of lb_slice.
        If values contains only one element, all members are set to it.

        """

        indices = [lb_slice.x_indices, lb_slice.y_indices, lb_slice.z_indices]
        N = [len(x) for x in indices]

        if N[0] * N[1] * N[2] == 0:
            raise AttributeError("Cannot set properties of an empty LBSlice")

        value = np.copy(value)
        attribute_shape = lb_slice.get_values(
            *np.zeros((3, 1), dtype=int), attribute).shape[3:]
        target_shape = (*N, *attribute_shape)

        # broadcast if only one element was provided
        if value.shape == attribute_shape:
            value = np.ones(target_shape) * value

        if value.shape != target_shape:
            raise ValueError(
                f"Input-dimensions of {attribute} array {value.shape} does not match slice dimensions {target_shape}.")

        lb_slice.set_values(*indices, attribute, value)

    def get_attribute(lb_slice, attribute):
        """
        Getter function that copies attribute from every member of
        lb_slice into an array (if possible).

        """

        indices = [lb_slice.x_indices, lb_slice.y_indices, lb_slice.z_indices]
        N = [len(x) for x in indices]

        if N[0] * N[1] * N[2] == 0:
            return np.empty(0, dtype=type(None))

        return lb_slice.get_values(*indices, attribute)

    for attribute_name in dir(LBFluidRoutines):
        if attribute_name in dir(LBSlice) or not isinstance(
                getattr(LBFluidRoutines, attribute_name), type(LBFluidRoutines.density)):
            continue

        # synthesize a new property
        new_property = property(
            functools.partial(get_attribute, attribute=attribute_name),
            functools.partial(set_attribute, attribute=attribute_name),
            doc=getattr(LBFluidRoutines, attribute_name).__doc__ or f'{attribute_name} for a slice')
        # attach the property to LBSlice
        setattr(LBSlice, attribute_name, new_property)


_add_lb_slice_properties()
