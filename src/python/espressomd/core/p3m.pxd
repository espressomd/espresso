include "myconfig.pxi"
IF ELECTROSTATICS:
    IF P3M:
        from core.p3m_common cimport P3MParameters
        cdef extern from "electrostatics_magnetostatics/p3m.hpp":
            int p3m_set_params(double r_cut, int * mesh, int cao, double alpha, double accuracy)
            void p3m_set_tune_params(double r_cut, int mesh[3], int cao, double alpha, double accuracy, int n_interpol)
            int p3m_set_mesh_offset(double x, double y, double z)
            int p3m_set_eps(double eps)
            int p3m_set_ninterpol(int n)
            int p3m_adaptive_tune(char ** log)

            ctypedef struct p3m_data_struct:
                P3MParameters params

            # links intern C-struct with python object
            cdef extern p3m_data_struct p3m
