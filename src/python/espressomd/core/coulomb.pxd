include "myconfig.pxi"
IF ELECTROSTATICS:
    cdef extern from "electrostatics_magnetostatics/coulomb.hpp":
        cdef enum CoulombMethod:
                    COULOMB_NONE, \
                    COULOMB_DH, \
                    COULOMB_P3M, \
                    COULOMB_MMM1D, \
                    COULOMB_MMM2D, \
                    COULOMB_ELC_P3M, \
                    COULOMB_RF, \
                    COULOMB_P3M_GPU, \
                    COULOMB_MMM1D_GPU, \
                    COULOMB_SCAFACOS

        ctypedef struct Coulomb_parameters:
            double prefactor
            CoulombMethod method

        cdef extern Coulomb_parameters coulomb

    cdef extern from "electrostatics_magnetostatics/coulomb.hpp" namespace "Coulomb":
        int set_prefactor(double prefactor)
        void deactivate_method()
