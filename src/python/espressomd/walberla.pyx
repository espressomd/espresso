
from .utils cimport make_Vector3d

cdef class WalberlaBlockForest:
    def __init__(self, box_size, double agrid, int ghost_layers = 1):
        mpi_init_walberla_blockforest(make_Vector3d(box_size), agrid, ghost_layers)
