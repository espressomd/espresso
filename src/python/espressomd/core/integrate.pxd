from libcpp cimport bool

cdef extern from "integrate.hpp":
    double time_step
    extern int integ_switch
    extern double sim_time
    extern double verlet_reuse
    extern double skin
    extern bool set_py_interrupt
