from .utils cimport Vector3i

cdef extern from "grid.hpp":
    cdef cppclass NodeGrid:
        Vector3i get_node_grid()
        void set_node_grid(const Vector3i &)

    NodeGrid node_grid
