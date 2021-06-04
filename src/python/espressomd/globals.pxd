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
from .utils cimport Vector3i

cdef extern from "grid.hpp":
    void mpi_set_node_grid(const Vector3i & node_grid)

cdef extern from "nonbonded_interactions/nonbonded_interaction_data.hpp":
    extern int max_seen_particle_type
    double maximal_cutoff_bonded()
    double maximal_cutoff_nonbonded()

cdef extern from "rattle.hpp":
    extern int n_rigidbonds
