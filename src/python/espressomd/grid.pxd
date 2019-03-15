from utils cimport Vector3i, Vector3d

cdef extern from "grid.hpp":
    Vector3d box_l
    Vector3d local_box_l
    Vector3i node_grid
    int periodic
    double min_box_l
