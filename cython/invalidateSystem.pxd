cdef extern from "config.h":
    pass

cdef extern from "global.h":
    cdef extern int INVALIDATE_SYSTEM

cdef extern from "communication.h":
    cdef void mpi_bcast_event(int event)
