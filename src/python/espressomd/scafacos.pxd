# Interface to the scafacos library. These are the methods shared between
# dipolar and electrostatics methods

from __future__ import print_function, absolute_import
include "myconfig.pxi"

from libcpp.string cimport string
from libcpp cimport bool
from libcpp.list cimport list
IF SCAFACOS:
    cdef extern from "electrostatics_magnetostatics/scafacos.hpp" namespace "Scafacos":
        void set_parameters(string & method_name, string & params, bool dipolar)
        string get_method_and_parameters()
        list[string] available_methods_core "Scafacos::available_methods" ()
        void free_handle()
