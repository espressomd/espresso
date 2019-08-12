from libcpp cimport bool

from utils cimport Vector3i, Vector3d

cdef extern from "grid.hpp":
    Vector3i node_grid

    cppclass BoxGeometry:
        void set_periodic(unsigned coord, bool value)
        bool periodic(unsigned coord)
        const Vector3d & length()
        void set_length(Vector3d)

    BoxGeometry box_geo

    Vector3d get_mi_vector(Vector3d, Vector3d, const BoxGeometry & )
    Vector3d folded_position(Vector3d, const BoxGeometry &)
    Vector3d unfolded_position(Vector3d, Vector3i, const Vector3d & )
