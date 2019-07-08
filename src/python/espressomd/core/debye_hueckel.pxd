include "myconfig.pxi"

IF ELECTROSTATICS:
    cdef extern from "electrostatics_magnetostatics/debye_hueckel.hpp":
        ctypedef struct Debye_hueckel_params:
            double r_cut
            double kappa

        cdef extern Debye_hueckel_params dh_params

        int dh_set_params(double kappa, double r_cut)
