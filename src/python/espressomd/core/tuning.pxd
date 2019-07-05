cdef extern from "tuning.hpp":
    cdef void c_tune_skin "tune_skin" (double min, double max, double tol, int steps)
