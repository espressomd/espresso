
include "myconfig.pxi"
from .utils cimport Vector3d

cdef extern from "grid_based_algorithms/walberla_blockforest.hpp":
    void mpi_init_walberla_blockforest(const Vector3d &box_size, double agrid, int n_ghost_layers) except +