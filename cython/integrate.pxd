cdef extern from "config.h":
    pass

cdef extern from "integrate.h":
    cdef void integrate_vv(int n_steps)
