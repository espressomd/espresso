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

from .utils cimport Vector3i

cdef extern from "grid_based_algorithms/lb_interface.hpp":

    double lb_lbfluid_get_tau() except +
    double lb_lbfluid_get_agrid() except +
    void lb_lbfluid_save_checkpoint(string filename, bool binary) except +
    void lb_lbfluid_load_checkpoint(string filename, bool binary) except +
    void lb_lbnode_remove_from_boundary(const Vector3i & ind) except +
    void lb_lbfluid_clear_boundaries() except +
    void lb_lbfluid_update_boundary_from_shape(const vector[int] & raster,
                                               const vector[double] & vel) except +
    void lb_lbfluid_update_boundary_from_list(const vector[int] & nodes_flat,
                                              const vector[double] & vel_flat) except +
    double lb_lbfluid_get_kT() except +
    double lb_lbfluid_get_lattice_speed() except +

cdef extern from "grid_based_algorithms/lb_particle_coupling.hpp":
    void lb_lbcoupling_set_rng_state(stdint.uint64_t) except +
    stdint.uint64_t lb_lbcoupling_get_rng_state() except +
    void lb_lbcoupling_set_gamma(double) except +
    double lb_lbcoupling_get_gamma() except +
    bool lb_lbcoupling_is_seed_required() except +
    void mpi_bcast_lb_particle_coupling()

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
