from libcpp.string cimport string
from libcpp.vector cimport vector

cdef extern from "random.hpp" namespace "Random":
    string mpi_random_get_stat()
    void mpi_random_set_stat(const vector[string] & stat)
    int get_state_size_of_generator()
