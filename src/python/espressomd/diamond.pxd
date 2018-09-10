from __future__ import print_function, absolute_import
include "myconfig.pxi"
from .c_analyze cimport PartCfg, partCfg

cdef extern from "polymer.hpp":
    int diamondC(PartCfg &, double a, double bond_length, int MPC, int N_CI, double val_nodes, double val_cM, double val_CI, int cM_dist, int nonet)


cdef class Diamond:
    cdef _params
