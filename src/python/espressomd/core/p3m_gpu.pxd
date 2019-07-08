include "myconfig.pxi"

IF ELECTROSTATICS:
    IF P3M:
        IF CUDA:
            cdef extern from "electrostatics_magnetostatics/p3m_gpu.hpp":
                void p3m_gpu_init(int cao, int * mesh, double alpha)

