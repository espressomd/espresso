from script_interface import  *
from script_interface cimport  *


cdef extern from "shapes/Shape.hpp" namespace "Shapes":
    cdef cppclass Shape:
        pass

cdef class PShape(PScriptInterface):
    cdef unique_ptr[Shape] si

cdef class PWall(PShape):
    pass
