from utils cimport Vector3i

cdef extern from "grid.hpp":
    double box_l[3]
    double local_box_l[3]
    extern Vector3i node_grid
    extern int periodic
    extern double min_box_l
