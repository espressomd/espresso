from __future__ import print_function, absolute_import
include "myconfig.pxi"
from .c_analyze cimport PartCfg, partCfg
from libcpp cimport bool

cdef extern from "diamond.hpp":
    int create_diamond(PartCfg &, double a, double bond_length, int MPC, int N_CI, double val_nodes, double val_cM, double val_CI, int cM_dist, bool nonet)


cdef class Diamond:
    cdef _params
