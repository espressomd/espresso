from libcpp.vector cimport vector
include "myconfig.pxi"

cdef extern from "communication.hpp":
    void mpi_bcast_cell_structure(int cs)
    int n_nodes
    vector[int] mpi_resort_particles(int global_flag)

IF ELECTROSTATICS and P3M:
    cdef extern from "communication.hpp":
        int mpi_iccp3m_init()
