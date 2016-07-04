include "myconfig.pxi"

IF P3M == 1 or DP3M == 1:
    cdef extern from "p3m-common.hpp":
        ctypedef struct p3m_parameter_struct:
            double alpha_L
            double r_cut_iL
            int    mesh[3]
            double mesh_off[3]
            int    cao
            int    inter
            double accuracy
            double epsilon
            double cao_cut[3]
            double a[3]
            double ai[3]
            double alpha
            double r_cut
            int    inter2
            int    cao3
            double additional_mesh[3]
