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

from libcpp cimport bool
from libcpp.vector cimport vector
from libcpp.string cimport string
from libc cimport stdint

from .actors cimport Actor
from .utils cimport Vector3d
from .utils cimport Vector3i
from .utils cimport Vector6d
from .utils cimport Vector19d

cdef class HydrodynamicInteraction(Actor):
    pass

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
    void lb_lbfluid_sanity_checks() except +
    void lb_lbfluid_print_vtk_velocity(string filename) except +
    void lb_lbfluid_print_vtk_velocity(string filename, vector[int] bb1, vector[int] bb2) except +
    void lb_lbfluid_print_vtk_boundary(string filename) except +
    void lb_lbfluid_print_velocity(string filename) except +
    void lb_lbfluid_print_boundary(string filename) except +
    void lb_lbfluid_save_checkpoint(string filename, bool binary) except +
    void lb_lbfluid_load_checkpoint(string filename, bool binary) except +
    void lb_lbfluid_set_lattice_switch(ActiveLB local_lattice_switch) except +
    Vector6d lb_lbfluid_get_stress() except +
    bool lb_lbnode_is_index_valid(const Vector3i & ind) except +
    Vector3i lb_lbfluid_get_shape() except +
    const Vector3d lb_lbnode_get_velocity(const Vector3i & ind) except +
    void lb_lbnode_set_velocity(const Vector3i & ind, const Vector3d & u) except +
    double lb_lbnode_get_density(const Vector3i & ind) except +
    void lb_lbnode_set_density(const Vector3i & ind, double density) except +
    const Vector6d lb_lbnode_get_stress(const Vector3i & ind) except +
    const Vector6d lb_lbnode_get_stress_neq(const Vector3i & ind) except +
    const Vector19d lb_lbnode_get_pop(const Vector3i & ind) except +
    void lb_lbnode_set_pop(const Vector3i & ind, const Vector19d & populations) except +
    int lb_lbnode_get_boundary(const Vector3i & ind) except +
    stdint.uint64_t lb_lbfluid_get_rng_state() except +
    void lb_lbfluid_set_rng_state(stdint.uint64_t) except +
    void lb_lbfluid_set_kT(double) except +
    double lb_lbfluid_get_kT() except +
    double lb_lbfluid_get_lattice_speed() except +
    void check_tau_time_step_consistency(double tau, double time_s) except +
    const Vector3d lb_lbfluid_get_interpolated_velocity(Vector3d & p) except +

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

###############################################
#
# Wrapper-functions for access to C-pointer: Set params
#
###############################################
cdef inline python_lbfluid_set_density(p_dens, agrid):
    cdef double c_dens

    # get pointers
    if isinstance(p_dens, float) or isinstance(p_dens, int):
        c_dens = <float > p_dens * agrid * agrid * agrid
    else:
        c_dens = p_dens * agrid * agrid * agrid
    # call c-function
    lb_lbfluid_set_density(c_dens)

###############################################

cdef inline python_lbfluid_set_viscosity(p_visc, p_agrid, p_tau):
    cdef double c_visc
    # get pointers
    if isinstance(p_visc, float) or isinstance(p_visc, int):
        c_visc = <float > p_visc * p_tau / (p_agrid * p_agrid)
    else:
        c_visc = p_visc * p_tau / (p_agrid * p_agrid)
    # call c-function
    lb_lbfluid_set_viscosity(c_visc)

###############################################

cdef inline python_lbfluid_set_agrid(p_agrid):
    cdef double c_agrid
    # get pointers
    c_agrid = p_agrid
    # call c-function
    lb_lbfluid_set_agrid(c_agrid)

###############################################

cdef inline python_lbfluid_set_bulk_viscosity(p_bvisc, p_agrid, p_tau):
    cdef double c_bvisc
    # get pointers
    if isinstance(p_bvisc, float) or isinstance(p_bvisc, int):
        c_bvisc = <float > p_bvisc * p_tau / (p_agrid * p_agrid)
    else:
        c_bvisc = p_bvisc * p_tau / (p_agrid * p_agrid)
    # call c-function
    lb_lbfluid_set_bulk_viscosity(c_bvisc)

###############################################

cdef inline python_lbfluid_set_gamma(p_gamma):
    cdef double c_gamma
    # get pointers
    if isinstance(p_gamma, float) or isinstance(p_gamma, int):
        c_gamma = <float > p_gamma
    else:
        c_gamma = p_gamma
    # call c-function
    lb_lbcoupling_set_gamma(c_gamma)

cdef inline python_lbfluid_set_gamma_odd(gamma_odd):
    cdef double c_gamma_odd
    # get pointers
    if isinstance(gamma_odd, float) or isinstance(gamma_odd, int):
        c_gamma_odd = <float > gamma_odd
    else:
        c_gamma_odd = gamma_odd
    # call c-function
    lb_lbfluid_set_gamma_odd(c_gamma_odd)

cdef inline python_lbfluid_set_gamma_even(gamma_even):
    cdef double c_gamma_even
    # get pointers
    if isinstance(gamma_even, float) or isinstance(gamma_even, int):
        c_gamma_even = <float > gamma_even
    else:
        c_gamma_even = gamma_even
    # call c-function
    lb_lbfluid_set_gamma_even(c_gamma_even)

###############################################
cdef inline python_lbfluid_set_ext_force_density(p_ext_force_density, p_agrid, p_tau):

    cdef Vector3d c_ext_force_density
    # unit conversion MD -> LB
    c_ext_force_density[
        0] = p_ext_force_density[
            0] * p_agrid * p_agrid * p_tau * p_tau
    c_ext_force_density[
        1] = p_ext_force_density[
            1] * p_agrid * p_agrid * p_tau * p_tau
    c_ext_force_density[
        2] = p_ext_force_density[
            2] * p_agrid * p_agrid * p_tau * p_tau
    lb_lbfluid_set_ext_force_density(c_ext_force_density)

cdef inline double python_lbfluid_get_density(agrid):
    return lb_lbfluid_get_density() / agrid / agrid / agrid

cdef inline double python_lbfluid_get_viscosity(p_agrid, p_tau):
    return lb_lbfluid_get_viscosity() / p_tau * (p_agrid * p_agrid)

cdef inline double python_lbfluid_get_bulk_viscosity(p_agrid, p_tau):
    return lb_lbfluid_get_bulk_viscosity() / p_tau * (p_agrid * p_agrid)

cdef inline double python_lbfluid_get_gamma():
    return lb_lbcoupling_get_gamma()

cdef inline Vector3d python_lbfluid_get_ext_force_density(p_agrid, p_tau):
    cdef Vector3d c_ext_force_density
    # call c-function
    c_ext_force_density = lb_lbfluid_get_ext_force_density()
    # unit conversion LB -> MD
    for i in range(3):
        c_ext_force_density[i] /= p_agrid * p_agrid * p_tau * p_tau
    return c_ext_force_density

cdef inline Vector6d python_lbfluid_get_stress(agrid, tau):
    cdef Vector6d stress = lb_lbfluid_get_stress()
    for i in range(6):
        stress[i] *= 1. / agrid * 1. / tau**2.0
    return stress

cdef inline void python_lbnode_set_velocity(Vector3i node, Vector3d velocity):
    cdef double inv_lattice_speed = lb_lbfluid_get_tau() / lb_lbfluid_get_agrid()
    cdef Vector3d c_velocity = velocity * inv_lattice_speed
    lb_lbnode_set_velocity(node, c_velocity)

cdef inline Vector3d python_lbnode_get_velocity(Vector3i node):
    cdef double lattice_speed = lb_lbfluid_get_agrid() / lb_lbfluid_get_tau()
    cdef Vector3d c_velocity = lb_lbnode_get_velocity(node)
    return c_velocity * lattice_speed

cdef inline void python_lbnode_set_density(Vector3i node, double density):
    cdef double agrid = lb_lbfluid_get_agrid()
    cdef double c_density = density * agrid * agrid * agrid
    lb_lbnode_set_density(node, c_density)

cdef inline double python_lbnode_get_density(Vector3i node):
    cdef double agrid = lb_lbfluid_get_agrid()
    cdef double c_density = lb_lbnode_get_density(node)
    return c_density / agrid / agrid / agrid

cdef inline Vector6d python_lbnode_get_stress(Vector3i node):
    cdef double tau = lb_lbfluid_get_tau()
    cdef double agrid = lb_lbfluid_get_agrid()
    cdef double unit_conversion = 1.0 / (tau * tau * agrid)
    cdef Vector6d c_stress = lb_lbnode_get_stress(node)
    return c_stress * unit_conversion

cdef inline Vector6d python_lbnode_get_stress_neq(Vector3i node):
    cdef double tau = lb_lbfluid_get_tau()
    cdef double agrid = lb_lbfluid_get_agrid()
    cdef double unit_conversion = 1.0 / (tau * tau * agrid)
    cdef Vector6d c_stress = lb_lbnode_get_stress_neq(node)
    return c_stress * unit_conversion
