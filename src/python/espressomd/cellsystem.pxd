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

from libcpp cimport bool
from libcpp.set cimport set as cpp_set
from libcpp.pair cimport pair, tuple
from libcpp.vector cimport vector
from .utils cimport Vector3i, Vector3d

cdef extern from "cells.hpp":
    cppclass PairInfo:
        int id1
        int id2
        Vector3d pos1
        Vector3d pos2
        Vector3d vec21
        int node

cdef extern from "communication.hpp":
    int n_nodes

cdef extern from "CellStructureType.hpp":
    ctypedef enum CellStructureType:
        CELL_STRUCTURE_REGULAR "CellStructureType::CELL_STRUCTURE_REGULAR"
        CELL_STRUCTURE_NSQUARE "CellStructureType::CELL_STRUCTURE_NSQUARE"
        CELL_STRUCTURE_HYBRID "CellStructureType::CELL_STRUCTURE_HYBRID"

cdef extern from "cells.hpp":
    ctypedef struct CellStructure:
        CellStructureType decomposition_type()
        bool use_verlet_list

    CellStructure cell_structure

    const RegularDecomposition * get_regular_decomposition()
    const HybridDecomposition * get_hybrid_decomposition()

    vector[pair[int, int]] mpi_get_pairs(double distance) except +
    vector[pair[int, int]] mpi_get_pairs_of_types(double distance, vector[int] types) except +
    vector[PairInfo] mpi_non_bonded_loop_trace()
    vector[int] mpi_resort_particles(int global_flag)
    void mpi_bcast_cell_structure(int cs)
    void mpi_set_hybrid_decomposition(cpp_set[int] n_square_types, double cutoff_regular)
    void mpi_set_use_verlet_lists(bool use_verlet_lists)

    pair[size_t, size_t] hybrid_parts_per_decomposition()

cdef extern from "tuning.hpp":
    cdef void c_tune_skin "tune_skin" (double min_skin, double max_skin, double tol, int int_steps, bool adjust_max_skin)

cdef extern from "integrate.hpp":
    extern double skin
    void mpi_set_skin(double skin)
    double get_verlet_reuse()

cdef extern from "RegularDecomposition.hpp":
    cppclass RegularDecomposition:
        Vector3i cell_grid
        double cell_size[3]

cdef extern from "HybridDecomposition.hpp":
    cppclass HybridDecomposition:
        Vector3i get_cell_grid()
        Vector3d get_cell_size()
        cpp_set[int] get_n_square_types()
        double get_cutoff_regular()

cdef extern from "grid.hpp":
    void mpi_set_node_grid(const Vector3i & node_grid)

cdef extern from "bonded_interactions/bonded_interaction_data.hpp":
    double maximal_cutoff_bonded()

cdef extern from "nonbonded_interactions/nonbonded_interaction_data.hpp":
    double maximal_cutoff_nonbonded()
