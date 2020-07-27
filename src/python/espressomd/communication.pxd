from libcpp.memory cimport shared_ptr
from boost cimport environment

cdef extern from "MpiCallbacks.hpp" namespace "Communication":
    cppclass MpiCallbacks:
        pass

cdef extern from "communication.hpp":
    shared_ptr[environment] mpi_init()
    void mpi_loop()
    int this_node

cdef extern from "communication.hpp" namespace "Communication":
    MpiCallbacks & mpiCallbacks()
    void init(shared_ptr[environment])
