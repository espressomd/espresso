include "myconfig.pxi"


IF DIPOLES==1:
    cdef extern from "interaction_data.hpp":
        ctypedef enum DipolarInteraction: 
            DIPOLAR_NONE = 0,
            DIPOLAR_P3M,
            DIPOLAR_MDLC_P3M,
            DIPOLAR_ALL_WITH_ALL_AND_NO_REPLICA,
            DIPOLAR_DS,
            DIPOLAR_MDLC_DS
        
        int dipolar_set_Dbjerrum(double bjerrum)
        
        ctypedef struct Coulomb_parameters:
            double Dprefactor
            DipolarInteraction Dmethod
        cdef extern Coulomb_parameters coulomb

        
    cdef extern from "magnetic_non_p3m_methods.hpp":
        int dawaanr_set_params()
        int mdds_set_params(int n_cut)
        int Ncut_off_magnetic_dipolar_direct_sum

