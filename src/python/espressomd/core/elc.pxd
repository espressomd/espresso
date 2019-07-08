from libcpp cimport bool
include "myconfig.pxi"

IF ELECTROSTATICS and P3M:
    cdef extern from "electrostatics_magnetostatics/elc.hpp":
        ctypedef struct ELC_struct:
            double maxPWerror
            double gap_size
            double far_cut
            int neutralize
            double delta_mid_top,
            double delta_mid_bot,
            bool const_pot,
            double pot_diff

        int ELC_set_params(double maxPWerror, double min_dist, double far_cut,
                           int neutralize, double delta_mid_top, double delta_mid_bot, bool const_pot, double pot_diff)

        # links intern C-struct with python object
        ELC_struct elc_params
