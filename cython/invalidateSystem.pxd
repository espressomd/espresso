cdef extern from "config.hpp":
    pass

cdef extern from "global.hpp":
    cdef extern int INVALIDATE_SYSTEM

cdef extern from "communication.hpp":
    cdef void mpi_bcast_event(int event)
