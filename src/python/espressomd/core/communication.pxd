from libcpp.vector cimport vector

cdef extern from "communication.hpp":
    void mpi_bcast_cell_structure(int cs)
    int n_nodes
    vector[int] mpi_resort_particles(int global_flag)
