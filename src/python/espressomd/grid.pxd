from libcpp cimport bool

from utils cimport Vector3i, Vector3d

cdef extern from "grid.hpp":
    Vector3d local_box_l
    Vector3i node_grid
    double min_box_l

    cppclass BoxGeometry:
        void set_periodic(unsigned coord, bool value)
        bool periodic(unsigned coord)
        const Vector3d & length()
        void set_length(Vector3d)

    BoxGeometry box_geo

    Vector3d get_mi_vector(const Vector3d & , const Vector3d & , const BoxGeometry & )
    Vector3d folded_position(const Vector3d &, const BoxGeometry & )
    Vector3d unfolded_position(const Vector3d & , const Vector3i & , const Vector3d & )
