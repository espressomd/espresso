include "myconfig.pxi"
from .analyze cimport PartCfg, partCfg

cdef extern from "diamond.hpp":
    int create_diamond(PartCfg &, double a, double bond_length, int MPC, int N_CI, double val_nodes, double val_cM, double val_CI, int cM_dist, int nonet)


cdef class Diamond:
    cdef _params
