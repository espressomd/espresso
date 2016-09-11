from __future__ import print_function, absolute_import
include "myconfig.pxi"

cdef extern from "polymer.hpp":
    int diamondC(double a, double bond_length, int MPC, int N_CI, double val_nodes, double val_cM, double val_CI, int cM_dist, int nonet)


cdef extern from "interaction_data.hpp":
    int n_bonded_ia

cdef class Diamond:
    cdef _params
