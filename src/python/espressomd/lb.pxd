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
from .utils cimport make_array_locked


cdef class LBFluidRoutines:
    cdef Vector3i node

##############################################
#
# extern functions and structs
#
##############################################

IF LB_WALBERLA:
    cdef extern from "grid_based_algorithms/lb_interface.hpp":

        cdef enum OutputVTK:
            pass

    cdef extern from "walberla_bridge/LBWalberlaBase.hpp" namespace 'OutputVTK':

        cdef OutputVTK output_vtk_density 'OutputVTK::density'
        cdef OutputVTK output_vtk_velocity_vector 'OutputVTK::velocity_vector'
        cdef OutputVTK output_vtk_pressure_tensor 'OutputVTK::pressure_tensor'

cdef extern from "grid_based_algorithms/lb_interface.hpp":

    void lb_lbfluid_set_tau(double c_tau) except +
    double lb_lbfluid_get_tau() except +
    double lb_lbfluid_get_viscosity() except +
    double lb_lbfluid_get_agrid() except +
    void lb_lbfluid_set_ext_force_density(const Vector3d forcedensity) except +
    const Vector3d lb_lbfluid_get_ext_force_density() except +
    double lb_lbfluid_get_bulk_viscosity() except +
    void lb_lbfluid_save_checkpoint(string filename, bool binary) except +
    void lb_lbfluid_load_checkpoint(string filename, bool binary) except +
    void lb_lbfluid_create_vtk(unsigned delta_N, unsigned initial_count, unsigned flag_observables, const string identifier, const string base_folder, const string execution_folder) except +
    void lb_lbfluid_switch_vtk(const string vtk_uid, int status) except +
    void lb_lbfluid_write_vtk(const string vtk_uid) except +
    Vector6d lb_lbfluid_get_pressure_tensor() except +
    bool lb_lbnode_is_index_valid(const Vector3i & ind) except +
    Vector3i lb_lbfluid_get_shape() except +
    const Vector3d lb_lbnode_get_velocity(const Vector3i & ind) except +
    const Vector3d lb_lbnode_get_velocity_at_boundary(const Vector3i & ind) except +
    const Vector3d lb_lbnode_get_boundary_force(const Vector3i & ind) except +
    const Vector3d lb_lbnode_get_last_applied_force(const Vector3i & ind) except +
    void lb_lbnode_set_velocity(const Vector3i & ind, const Vector3d & u) except +
    void lb_lbnode_set_velocity_at_boundary(const Vector3i & ind, const Vector3d & u) except +
    void lb_lbnode_set_last_applied_force(const Vector3i & ind, const Vector3d & f) except +
    double lb_lbnode_get_density(const Vector3i & ind) except +
    void lb_lbnode_set_density(const Vector3i & ind, double density) except +
    const Vector6d lb_lbnode_get_pressure_tensor(const Vector3i & ind) except +
    const vector[double] lb_lbnode_get_pop(const Vector3i & ind) except +
    void lb_lbnode_set_pop(const Vector3i & ind, const vector[double] & populations) except +
    bool lb_lbnode_is_boundary(const Vector3i & ind) except +
    void lb_lbnode_remove_from_boundary(const Vector3i & ind) except +
    void lb_lbfluid_clear_boundaries() except +
    void lb_lbfluid_update_boundary_from_shape(const vector[int] & raster,
                                               const vector[double] & vel) except +
    void lb_lbfluid_update_boundary_from_list(const vector[int] & nodes_flat,
                                              const vector[double] & vel_flat) except +
    stdint.uint64_t lb_lbfluid_get_rng_state() except +
    void lb_lbfluid_set_rng_state(stdint.uint64_t) except +
    void lb_lbfluid_set_kT(double) except +
    double lb_lbfluid_get_kT() except +
    double lb_lbfluid_get_lattice_speed() except +
    void check_tau_time_step_consistency(double tau, double time_s) except +
    const Vector3d lb_lbfluid_get_interpolated_velocity(const Vector3d & p) except +
    void lb_lbfluid_add_force_at_pos(const Vector3d & p, const Vector3d & f) except +

cdef extern from "grid_based_algorithms/lb_particle_coupling.hpp":
    void lb_lbcoupling_set_rng_state(stdint.uint64_t) except +
    stdint.uint64_t lb_lbcoupling_get_rng_state() except +
    void lb_lbcoupling_set_gamma(double) except +
    double lb_lbcoupling_get_gamma() except +
    bool lb_lbcoupling_is_seed_required() except +
    void mpi_bcast_lb_particle_coupling()


##############################################
#
# Wrapper-functions to handle unit conversions
#
##############################################

cdef inline python_lbfluid_set_ext_force_density(Vector3d ext_force_density, double agrid, double tau) except +:
    lb_lbfluid_set_ext_force_density(ext_force_density * agrid**2 * tau**2)

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

cdef inline python_lbnode_set_velocity_at_boundary(Vector3i node, Vector3d velocity) except +:
    cdef double inv_lattice_speed = 1.0 / lb_lbfluid_get_lattice_speed()
    cdef Vector3d c_velocity = velocity * inv_lattice_speed
    lb_lbnode_set_velocity_at_boundary(node, c_velocity)

cdef inline python_lbnode_get_velocity(Vector3i node) except +:
    cdef Vector3d c_velocity = lb_lbnode_get_velocity(node)
    return make_array_locked(c_velocity * lb_lbfluid_get_lattice_speed())

cdef inline python_lbnode_get_interpolated_velocity(Vector3d pos) except +:
    cdef Vector3d c_velocity = lb_lbfluid_get_interpolated_velocity(pos)
    return make_array_locked(c_velocity * lb_lbfluid_get_lattice_speed())

cdef inline python_lbnode_set_last_applied_force(Vector3i node, Vector3d force) except +:
    cdef double unit_conversion = lb_lbfluid_get_tau()**2 / lb_lbfluid_get_agrid()
    cdef Vector3d c_f = force * unit_conversion
    lb_lbnode_set_last_applied_force(node, c_f)

cdef inline python_lbnode_get_velocity_at_boundary(Vector3i node) except +:
    cdef Vector3d c_velocity = lb_lbnode_get_velocity_at_boundary(node)
    cdef double lattice_speed = lb_lbfluid_get_lattice_speed()
    return make_array_locked(c_velocity * lattice_speed)

cdef inline python_lbnode_get_boundary_force(Vector3i node) except +:
    cdef Vector3d force = lb_lbnode_get_boundary_force(node)
    cdef double unit_conversion = lb_lbfluid_get_agrid() / lb_lbfluid_get_tau()**2
    return make_array_locked(force * unit_conversion)

cdef inline python_lbnode_get_last_applied_force(Vector3i node) except +:
    cdef Vector3d c_f = lb_lbnode_get_last_applied_force(node)
    cdef double unit_conversion = lb_lbfluid_get_agrid() / lb_lbfluid_get_tau()**2
    return make_array_locked(c_f * unit_conversion)

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

cdef inline python_lb_lbfluid_update_boundary_from_shape(int[:] raster_view, double[:] vel_view) except +:
    cdef vector[int] raster
    cdef vector[double] vel
    cdef int * raster_ptr = &raster_view[0]
    cdef double * vel_ptr = &vel_view[0]
    raster.assign(raster_ptr, raster_ptr + len(raster_view))
    vel.assign(vel_ptr, vel_ptr + len(vel_view))
    lb_lbfluid_update_boundary_from_shape(raster, vel)

cdef inline python_lb_lbfluid_update_boundary_from_list(int[:] nodes_view, double[:] vel_view) except +:
    cdef vector[int] nodes
    cdef vector[double] vel
    cdef int * nodes_ptr = &nodes_view[0]
    cdef double * vel_ptr = &vel_view[0]
    nodes.assign(nodes_ptr, nodes_ptr + len(nodes_view))
    vel.assign(vel_ptr, vel_ptr + len(vel_view))
    lb_lbfluid_update_boundary_from_list(nodes, vel)
