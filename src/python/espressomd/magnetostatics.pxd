from __future__ import print_function, absolute_import
include "myconfig.pxi"


IF DIPOLES == 1:
    cdef extern from "communication.hpp":
        void mpi_bcast_coulomb_params()

    cdef extern from "nonbonded_interactions/nonbonded_interaction_data.hpp":
        ctypedef enum dipolar_interaction "DipolarInteraction":
            DIPOLAR_NONE = 0,
            DIPOLAR_P3M,
            DIPOLAR_MDLC_P3M,
            DIPOLAR_DS,
            DIPOLAR_MDLC_DS,
            DIPOLAR_SCAFACOS

        int dipolar_set_Dprefactor(double prefactor)

        ctypedef struct coulomb_parameters "Coulomb_parameters":
            double Dprefactor
            dipolar_interaction Dmethod
        cdef extern coulomb_parameters coulomb

    cdef extern from "electrostatics_magnetostatics/mdds.hpp":
        void mdds_set_params(int n_cut)
        int mdds_n_replicas

    IF(CUDA == 1) and (ROTATION == 1):
        cdef extern from "actor/DipolarDirectSum.hpp":
            void activate_dipolar_direct_sum_gpu()
            void deactivate_dipolar_direct_sum_gpu()

    IF(DIPOLAR_BARNES_HUT == 1):
        cdef extern from "actor/DipolarBarnesHut.hpp":
            void activate_dipolar_barnes_hut(float epssq, float itolsq)
            void deactivate_dipolar_barnes_hut()

IF DP3M == 1:
    from p3m_common cimport p3m_parameter_struct

    cdef extern from "electrostatics_magnetostatics/p3m-dipolar.hpp":
        int dp3m_set_params(double r_cut, int mesh, int cao, double alpha, double accuracy)
        void dp3m_set_tune_params(double r_cut, int mesh, int cao, double alpha, double accuracy, int n_interpol)
        int dp3m_set_mesh_offset(double x, double y, double z)
        int dp3m_set_eps(double eps)
        int dp3m_set_ninterpol(int n)
        int dp3m_adaptive_tune(char ** log)
        int dp3m_deactivate()

        ctypedef struct dp3m_data_struct:
            p3m_parameter_struct params

        # links intern C-struct with python object
        cdef extern dp3m_data_struct dp3m

        # Convert C arguments into numpy array
