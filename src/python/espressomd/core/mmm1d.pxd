include "myconfig.pxi"

IF ELECTROSTATICS:
    cdef extern from "electrostatics_magnetostatics/mmm1d.hpp":
        ctypedef struct MMM1D_struct:
            double far_switch_radius_2;
            double maxPWerror;
            int    bessel_cutoff;

        cdef extern MMM1D_struct mmm1d_params;

        int MMM1D_set_params(double switch_rad, double maxPWerror);
        void MMM1D_init();
        int MMM1D_sanity_checks();
        int mmm1d_tune(char ** log);


