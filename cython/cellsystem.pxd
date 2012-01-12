cdef extern from "communication.h":
    void mpi_bcast_parameter(int p)
