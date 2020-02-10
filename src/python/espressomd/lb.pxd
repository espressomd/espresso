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

    #
    #
    # extern functions and structs
    #
    #

cdef extern from "grid_based_algorithms/lb_interface.hpp" namespace "ActiveLB":
    cdef ActiveLB NONE
    cdef ActiveLB WALBERLA 

cdef extern from "grid_based_algorithms/lb_interface.hpp":

    cdef enum ActiveLB:
        pass
    void lb_lbfluid_set_tau(double c_tau) except +
    double lb_lbfluid_get_tau() except +
    double lb_lbfluid_get_viscosity() except +
    double lb_lbfluid_get_agrid() except +
    void lb_lbfluid_set_ext_force_density(const Vector3d forcedensity) except +
    const Vector3d lb_lbfluid_get_ext_force_density() except +
    double lb_lbfluid_get_bulk_viscosity() except +
    void lb_lbfluid_sanity_checks() except +
    void lb_lbfluid_print_vtk_velocity(string filename) except +
    void lb_lbfluid_print_vtk_velocity(string filename, vector[int] bb1, vector[int] bb2) except +
    void lb_lbfluid_print_vtk_boundary(string filename) except +
    void lb_lbfluid_print_velocity(string filename) except +
    void lb_lbfluid_print_boundary(string filename) except +
    void lb_lbfluid_save_checkpoint(string filename, int binary) except +
    void lb_lbfluid_load_checkpoint(string filename, int binary) except +
    void lb_lbfluid_set_lattice_switch(ActiveLB local_lattice_switch) except +
    ActiveLB lb_lbfluid_get_lattice_switch() except +
    Vector6d lb_lbfluid_get_stress() except +
    bool lb_lbnode_is_index_valid(const Vector3i & ind) except +
    Vector3i lb_lbfluid_get_shape() except +
    const Vector3d lb_lbnode_get_velocity(const Vector3i & ind) except +
    void lb_lbnode_set_velocity(const Vector3i & ind, const Vector3d & u) except +
    double lb_lbnode_get_density(const Vector3i & ind) except +
    void lb_lbnode_set_density(const Vector3i & ind, double density) except +
    const Vector6d lb_lbnode_get_stress(const Vector3i & ind) except +
    const Vector19d lb_lbnode_get_pop(const Vector3i & ind) except +
    void lb_lbnode_set_pop(const Vector3i & ind, const Vector19d & populations) except +
    bool lb_lbnode_is_boundary(const Vector3i & ind) except +
    stdint.uint64_t lb_lbfluid_get_rng_state() except +
    void lb_lbfluid_set_rng_state(stdint.uint64_t) except +
    void lb_lbfluid_set_kT(double) except +
    double lb_lbfluid_get_kT() except +
    double lb_lbfluid_get_lattice_speed() except +
    void check_tau_time_step_consistency(double tau, double time_s) except +
    const Vector3d lb_lbfluid_get_interpolated_velocity(Vector3d & p) except +

cdef extern from "grid_based_algorithms/lb_particle_coupling.hpp":
    void lb_lbcoupling_set_rng_state(stdint.uint64_t) except +
    stdint.uint64_t lb_lbcoupling_get_rng_state() except +
    void lb_lbcoupling_set_gamma(double) except +
    double lb_lbcoupling_get_gamma() except +
    bool lb_lbcoupling_is_seed_required() except +

cdef extern from "grid_based_algorithms/lb_interpolation.hpp":
    cdef cppclass InterpolationOrder:
        pass
    void lb_lbinterpolation_set_interpolation_order(InterpolationOrder & order) except +

cdef extern from "grid_based_algorithms/lb_interpolation.hpp" namespace "InterpolationOrder":
    cdef InterpolationOrder linear
    cdef InterpolationOrder quadratic

IF LB_WALBERLA:
    cdef extern from "grid_based_algorithms/lb_walberla_instance.hpp":
        void mpi_init_lb_walberla(double viscosity, double density, double agrid, double tau) except +
        void mpi_destruct_lb_walberla() except +


#
#
# Wrapper-functions for access to C-pointer: Set params
#
#

cdef inline python_lbfluid_set_ext_force_density(p_ext_force_density, p_agrid, p_tau) except +:

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

cdef inline double python_lbfluid_get_viscosity(p_agrid, p_tau) except +:
    return lb_lbfluid_get_viscosity() / p_tau * (p_agrid * p_agrid)

cdef inline double python_lbfluid_get_bulk_viscosity(p_agrid, p_tau) except +:
    return lb_lbfluid_get_bulk_viscosity() / p_tau * (p_agrid * p_agrid)

cdef inline double python_lbfluid_get_gamma() except +:
    return lb_lbcoupling_get_gamma()

cdef inline Vector3d python_lbfluid_get_ext_force_density(p_agrid, p_tau) except +:
    cdef Vector3d c_ext_force_density
    # call c-function
    c_ext_force_density = lb_lbfluid_get_ext_force_density()
    # unit conversion LB -> MD
    for i in range(3):
        c_ext_force_density[i] /= p_agrid * p_agrid * p_tau * p_tau
    return c_ext_force_density

cdef inline Vector6d python_lbfluid_get_stress(agrid, tau) except +:
    cdef Vector6d stress = lb_lbfluid_get_stress()
    for i in range(6):
        stress[i] *= 1. / agrid * 1. / tau**2.0
    return stress

cdef inline void python_lbnode_set_velocity(Vector3i node, Vector3d velocity) except +:
    cdef double inv_lattice_speed = lb_lbfluid_get_tau() / lb_lbfluid_get_agrid()
    cdef Vector3d c_velocity = velocity * inv_lattice_speed
    lb_lbnode_set_velocity(node, c_velocity)

cdef inline Vector3d python_lbnode_get_velocity(Vector3i node) except +:
    cdef double lattice_speed = lb_lbfluid_get_agrid() / lb_lbfluid_get_tau()
    cdef Vector3d c_velocity = lb_lbnode_get_velocity(node)
    return c_velocity * lattice_speed

cdef inline Vector6d python_lbnode_get_stress(Vector3i node) except +:
    cdef double tau = lb_lbfluid_get_tau()
    cdef double agrid = lb_lbfluid_get_agrid()
    cdef double unit_conversion = 1.0 / (tau * tau * agrid)
    cdef Vector6d c_stress = lb_lbnode_get_stress(node)
    return c_stress * unit_conversion

cdef inline double python_lbnode_get_density(Vector3i node) except +: 
    cdef double agrid = lb_lbfluid_get_agrid()
    cdef double c_density = lb_lbnode_get_density(node)
    return c_density / agrid / agrid / agrid 

cdef inline void python_lbnode_set_density(Vector3i node, double density) except +:
    cdef double agrid = lb_lbfluid_get_agrid()
    cdef double c_density = density * agrid * agrid * agrid
    lb_lbnode_set_density(node, c_density)

cdef inline python_lbfluid_set_gamma(p_gamma) except +:
    cdef double c_gamma
    # get pointers
    if isinstance(p_gamma, float) or isinstance(p_gamma, int):
        c_gamma = <float > p_gamma
    else:
        c_gamma = p_gamma
        # call c-function
        lb_lbcoupling_set_gamma(c_gamma)
