# Interface to the scafacos library. These are the methods shared between
# dipolar and electrostatics methods

from __future__ import print_function, absolute_import
include "myconfig.pxi"

from libcpp.string cimport string
from libcpp cimport bool
from libcpp.list cimport list
IF SCAFACOS == 1:
    cdef extern from "scafacos.hpp" namespace "Scafacos":
        cdef void set_parameters(string & method_name, string & params, bool dipolar)
        cdef string get_parameters()
        cpdef list[string] available_methods()
