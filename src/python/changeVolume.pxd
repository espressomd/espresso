cdef extern from "config.hpp":
    pass

cdef extern from "grid.hpp":
    cdef void rescale_boxl(int dir, double d_new)

