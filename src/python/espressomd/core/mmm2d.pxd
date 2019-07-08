include "myconfig.pxi"
from libcpp cimport bool

IF ELECTROSTATICS:
    cdef extern from "electrostatics_magnetostatics/mmm2d.hpp":
        ctypedef struct MMM2D_struct:
            double maxPWerror;
            double far_cut;
            double far_cut2;
            int far_calculated;
            bool dielectric_contrast_on;
            bool const_pot_on;
            double pot_diff;
            double delta_mid_top;
            double delta_mid_bot;
            double delta_mult;

        cdef extern MMM2D_struct mmm2d_params;

        int MMM2D_set_params(double maxPWerror, double far_cut, double delta_top, double delta_bot, bool const_pot_on, double pot_diff);

        void MMM2D_init();

        int MMM2D_sanity_checks();
