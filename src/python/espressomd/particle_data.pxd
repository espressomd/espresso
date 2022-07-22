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
from libcpp cimport bool
from libcpp.vector cimport vector  # import std::vector as vector

include "myconfig.pxi"

cdef extern from "particle_node.hpp":
    void prefetch_particle_data(vector[int] ids)

    bool particle_exists(int p_id)

    int get_maximal_particle_id()

cdef class _ParticleSliceImpl:
    cdef public id_selection
    cdef int _chunk_size
