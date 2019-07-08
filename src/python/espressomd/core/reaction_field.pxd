include "myconfig.pxi"

IF ELECTROSTATICS:
     cdef extern from "electrostatics_magnetostatics/reaction_field.hpp":
        ctypedef struct Reaction_field_params:
            double kappa
            double epsilon1
            double epsilon2
            double r_cut

        cdef extern Reaction_field_params rf_params

        int rf_set_params(double kappa, double epsilon1, double epsilon2,
                          double r_cut)

