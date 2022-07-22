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
# Here we create something to handle particles
from .utils cimport Vector4d, Vector3d, Vector3i
from libcpp cimport bool
from libcpp.vector cimport vector  # import std::vector as vector

include "myconfig.pxi"

cdef extern from "particle_node.hpp":
    void prefetch_particle_data(vector[int] ids)

    void place_particle(int p_id, const Vector3d & pos) except +

    void remove_particle(int p_id) except +

    void remove_all_particles() except +

    bool particle_exists(int p_id)

    int get_particle_node(int p_id) except +

    vector[int] get_particle_ids() except +

    int get_maximal_particle_id()
    int get_n_part()

cdef extern from "nonbonded_interactions/nonbonded_interaction_data.hpp":
    int max_seen_particle_type

cdef extern from "bonded_interactions/rigid_bond.hpp":
    extern int n_rigidbonds

cdef class _ParticleSliceImpl:
    cdef public id_selection
    cdef int _chunk_size
