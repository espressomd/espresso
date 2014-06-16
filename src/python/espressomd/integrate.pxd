cdef extern from "config.hpp":
    pass

cdef extern from "integrate.hpp":
    cdef void integrate_vv(int n_steps, int reuse_forces)
