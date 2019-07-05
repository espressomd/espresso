from core.utils cimport Vector9d

cdef extern from "dpd.hpp":
    Vector9d dpd_stress()
