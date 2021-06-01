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
from .utils cimport Vector3d, Vector3i

cdef extern from "grid.hpp":
    void mpi_set_box_length(Vector3d length) except +
    void mpi_set_periodicity(bool x, bool y, bool z)
    void mpi_set_node_grid(const Vector3i & node_grid)

cdef extern from "integrate.hpp":
    double get_time_step()
    extern int integ_switch
    double get_sim_time()
    extern double verlet_reuse
    extern double skin
    void mpi_set_time_step(double time_step) except +
    void mpi_set_skin(double skin)
    void mpi_set_time(double time)

cdef extern from "nonbonded_interactions/nonbonded_interaction_data.hpp":
    extern int max_seen_particle_type
    extern double min_global_cut
    double maximal_cutoff_bonded()
    double maximal_cutoff_nonbonded()
    void mpi_set_min_global_cut(double min_global_cut)

cdef extern from "rattle.hpp":
    extern int n_rigidbonds

cdef extern from "tuning.hpp":
    extern int timing_samples

cdef extern from "object-in-fluid/oif_global_forces.hpp":
    int max_oif_objects
    void mpi_set_max_oif_objects(int max_oif_objects)

cdef extern from "forcecap.hpp":
    double forcecap_get()
    void mpi_set_forcecap(double forcecap)
