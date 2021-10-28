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

from .utils cimport Vector3d

cdef extern from "grid.hpp":
    void mpi_set_box_length(Vector3d length) except +
    void mpi_set_periodicity(bool x, bool y, bool z)
    void rescale_boxl(int dir, double d_new)

cdef extern from "rotate_system.hpp":
    void mpi_rotate_system(double phi, double theta, double alpha)

IF EXCLUSIONS:
    cdef extern from "particle_data.hpp":
        void auto_exclusions(int distance)

cdef extern from "particle_data.hpp":
    int init_type_map(int type) except +
    int number_of_particles_with_type(int type) except +

cdef extern from "object-in-fluid/oif_global_forces.hpp":
    int max_oif_objects
    void mpi_set_max_oif_objects(int max_oif_objects)

cdef extern from "nonbonded_interactions/nonbonded_interaction_data.hpp":
    extern double min_global_cut
    void mpi_set_min_global_cut(double min_global_cut)
