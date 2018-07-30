
cdef extern from "tuning.hpp":
    cdef void tune_skin(double min, double max, double tol, int steps)
