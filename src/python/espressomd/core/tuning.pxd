cdef extern from "tuning.hpp":
    extern int timing_samples
    cdef void tune_skin(double min, double max, double tol, int steps)
