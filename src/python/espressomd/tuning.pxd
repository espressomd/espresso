from libcpp cimport bool

cdef extern from "tuning.hpp":
    cdef void tune_skin(double min_skin, double max_skin, double tol, int int_steps, bool adjust_max_skin)
