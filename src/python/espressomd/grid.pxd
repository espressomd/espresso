from .utils cimport Vector3i

cdef extern from "grid.hpp":
    cdef cppclass NodeGrid:
        Vector3i get_node_grid()
        void set_node_grid(const Vector3i &)

    void mpi_bcast_node_grid()

    NodeGrid node_grid
