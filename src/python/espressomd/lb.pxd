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

from libcpp cimport bool
from libcpp.vector cimport vector
from libcpp.string cimport string
from libc cimport stdint

from .utils cimport Vector3d
from .utils cimport Vector3i
from .utils cimport Vector6d
from .utils cimport Vector19d
from .utils cimport make_array_locked

cdef class FluidActor:
    cdef public _isactive
    cdef public _params
    cdef public system

cdef class HydrodynamicInteraction(FluidActor):
    pass

cdef class LBFluidRoutines:
    cdef Vector3i node

##############################################
#
# extern functions and structs
#
##############################################

cdef extern from "grid_based_algorithms/lb_interface.hpp" namespace "ActiveLB":
    cdef ActiveLB NONE
    cdef ActiveLB CPU
    cdef ActiveLB GPU

cdef extern from "grid_based_algorithms/lb_interface.hpp":

    cdef enum ActiveLB:
        pass
    void lb_lbfluid_set_tau(double c_tau) except +
    double lb_lbfluid_get_tau() except +
    void lb_lbfluid_set_density(double c_dens) except +
    double lb_lbfluid_get_density() except +
    void lb_lbfluid_set_viscosity(double c_visc) except +
    double lb_lbfluid_get_viscosity() except +
    void lb_lbfluid_set_agrid(double c_agrid) except +
    double lb_lbfluid_get_agrid() except +
    void lb_lbfluid_set_gamma_odd(double c_gamma_odd) except +
    double lb_lbfluid_get_gamma_odd() except +
    void lb_lbfluid_set_gamma_even(double c_gamma_even) except +
    double lb_lbfluid_get_gamma_even() except +
    void lb_lbfluid_set_ext_force_density(const Vector3d forcedensity) except +
    const Vector3d lb_lbfluid_get_ext_force_density() except +
    void lb_lbfluid_set_bulk_viscosity(double c_bulk_visc) except +
    double lb_lbfluid_get_bulk_viscosity() except +
    void lb_lbfluid_print_vtk_velocity(string filename) except +
    void lb_lbfluid_print_vtk_velocity(string filename, vector[int] bb1, vector[int] bb2) except +
    void lb_lbfluid_print_vtk_boundary(string filename) except +
    void lb_lbfluid_print_velocity(string filename) except +
    void lb_lbfluid_print_boundary(string filename) except +
    void lb_lbfluid_save_checkpoint(string filename, bool binary) except +
    void lb_lbfluid_load_checkpoint(string filename, bool binary) except +
    void lb_lbfluid_set_lattice_switch(ActiveLB local_lattice_switch) except +
    Vector6d lb_lbfluid_get_pressure_tensor() except +
    bool lb_lbnode_is_index_valid(const Vector3i & ind) except +
    Vector3i lb_lbfluid_get_shape() except +
    const Vector3d lb_lbnode_get_velocity(const Vector3i & ind) except +
    void lb_lbnode_set_velocity(const Vector3i & ind, const Vector3d & u) except +
    double lb_lbnode_get_density(const Vector3i & ind) except +
    void lb_lbnode_set_density(const Vector3i & ind, double density) except +
    const Vector6d lb_lbnode_get_pressure_tensor(const Vector3i & ind) except +
    const Vector6d lb_lbnode_get_pressure_tensor_neq(const Vector3i & ind) except +
    const Vector19d lb_lbnode_get_pop(const Vector3i & ind) except +
    void lb_lbnode_set_pop(const Vector3i & ind, const Vector19d & populations) except +
    int lb_lbnode_get_boundary(const Vector3i & ind) except +
    stdint.uint64_t lb_lbfluid_get_rng_state() except +
    void lb_lbfluid_set_rng_state(stdint.uint64_t) except +
    void lb_lbfluid_set_kT(double) except +
    double lb_lbfluid_get_kT() except +
    double lb_lbfluid_get_lattice_speed() except +
    void check_tau_time_step_consistency(double tau, double time_s) except +
    const Vector3d lb_lbfluid_get_interpolated_velocity(const Vector3d & p) except +

cdef extern from "grid_based_algorithms/lb_particle_coupling.hpp":
    void lb_lbcoupling_set_rng_state(stdint.uint64_t)
    stdint.uint64_t lb_lbcoupling_get_rng_state() except +
    void lb_lbcoupling_set_gamma(double)
    double lb_lbcoupling_get_gamma() except +
    bool lb_lbcoupling_is_seed_required()

cdef extern from "grid_based_algorithms/lbgpu.hpp":
    void linear_velocity_interpolation(double * positions, double * velocities, int length)
    void quadratic_velocity_interpolation(double * positions, double * velocities, int length)

cdef extern from "grid_based_algorithms/lb_interpolation.hpp":
    cdef cppclass InterpolationOrder:
        pass
    void lb_lbinterpolation_set_interpolation_order(InterpolationOrder & order)

cdef extern from "grid_based_algorithms/lb_interpolation.hpp" namespace "InterpolationOrder":
    cdef InterpolationOrder linear
    cdef InterpolationOrder quadratic

cdef extern from "integrate.hpp":
    double get_time_step()

##############################################
#
# Wrapper-functions to handle unit conversions
#
##############################################

cdef inline python_lbfluid_set_density(double dens, double agrid) except +:
    lb_lbfluid_set_density(dens * agrid**3)

cdef inline python_lbfluid_set_viscosity(double visc, double agrid, double tau) except +:
    lb_lbfluid_set_viscosity(visc * tau / agrid**2)

cdef inline python_lbfluid_set_agrid(double agrid) except +:
    lb_lbfluid_set_agrid(agrid)

cdef inline python_lbfluid_set_bulk_viscosity(double bvisc, double agrid, double tau) except +:
    lb_lbfluid_set_bulk_viscosity(bvisc * tau / agrid**2)

cdef inline python_lbfluid_set_gamma(double gamma) except +:
    lb_lbcoupling_set_gamma(gamma)

cdef inline python_lbfluid_set_gamma_odd(double gamma_odd) except +:
    lb_lbfluid_set_gamma_odd(gamma_odd)

cdef inline python_lbfluid_set_gamma_even(double gamma_even) except +:
    lb_lbfluid_set_gamma_even(gamma_even)

cdef inline python_lbfluid_set_ext_force_density(Vector3d ext_force_density, double agrid, double tau) except +:
    lb_lbfluid_set_ext_force_density(ext_force_density * agrid**2 * tau**2)

cdef inline python_lbfluid_get_density(double agrid) except +:
    return lb_lbfluid_get_density() / agrid**3

cdef inline python_lbfluid_get_viscosity(double agrid, double tau) except +:
    return lb_lbfluid_get_viscosity() / tau * agrid**2

cdef inline python_lbfluid_get_bulk_viscosity(double agrid, double tau) except +:
    return lb_lbfluid_get_bulk_viscosity() / tau * agrid**2

cdef inline python_lbfluid_get_gamma() except +:
    return lb_lbcoupling_get_gamma()

cdef inline python_lbfluid_get_ext_force_density(double agrid, double tau) except +:
    cdef Vector3d ext_force_density = lb_lbfluid_get_ext_force_density()
    return make_array_locked(ext_force_density / (agrid**2 * tau**2))

cdef inline python_lbfluid_get_pressure_tensor(double agrid, double tau) except +:
    cdef Vector6d c_tensor = lb_lbfluid_get_pressure_tensor()
    cdef double unit_conversion = 1.0 / (agrid * tau**2)
    cdef Vector6d p_tensor = c_tensor * unit_conversion
    return [[p_tensor[0], p_tensor[1], p_tensor[3]],
            [p_tensor[1], p_tensor[2], p_tensor[4]],
            [p_tensor[3], p_tensor[4], p_tensor[5]]]

cdef inline python_lbnode_set_velocity(Vector3i node, Vector3d velocity) except +:
    lb_lbnode_set_velocity(node, velocity / lb_lbfluid_get_lattice_speed())

cdef inline python_lbnode_get_velocity(Vector3i node) except +:
    cdef Vector3d c_velocity = lb_lbnode_get_velocity(node)
    return make_array_locked(c_velocity * lb_lbfluid_get_lattice_speed())

cdef inline python_lbnode_get_interpolated_velocity(Vector3d pos) except +:
    cdef Vector3d c_velocity = lb_lbfluid_get_interpolated_velocity(pos)
    return make_array_locked(c_velocity * lb_lbfluid_get_lattice_speed())

cdef inline python_lbnode_set_density(Vector3i node, double density) except +:
    cdef double agrid = lb_lbfluid_get_agrid()
    lb_lbnode_set_density(node, density * agrid**3)

cdef inline python_lbnode_get_density(Vector3i node) except +:
    cdef double c_density = lb_lbnode_get_density(node)
    cdef double agrid = lb_lbfluid_get_agrid()
    return c_density / agrid**3

cdef inline python_lbnode_get_pressure_tensor(Vector3i node) except +:
    cdef Vector6d c_tensor = lb_lbnode_get_pressure_tensor(node)
    cdef double tau = lb_lbfluid_get_tau()
    cdef double agrid = lb_lbfluid_get_agrid()
    cdef double unit_conversion = 1.0 / (tau**2 * agrid)
    cdef Vector6d p_tensor = c_tensor * unit_conversion
    return [[p_tensor[0], p_tensor[1], p_tensor[3]],
            [p_tensor[1], p_tensor[2], p_tensor[4]],
            [p_tensor[3], p_tensor[4], p_tensor[5]]]

cdef inline python_lbnode_get_pressure_tensor_neq(Vector3i node) except +:
    cdef Vector6d c_tensor = lb_lbnode_get_pressure_tensor_neq(node)
    cdef double tau = lb_lbfluid_get_tau()
    cdef double agrid = lb_lbfluid_get_agrid()
    cdef double unit_conversion = 1.0 / (tau**2 * agrid)
    cdef Vector6d p_tensor = c_tensor * unit_conversion
    return [[p_tensor[0], p_tensor[1], p_tensor[3]],
            [p_tensor[1], p_tensor[2], p_tensor[4]],
            [p_tensor[3], p_tensor[4], p_tensor[5]]]
