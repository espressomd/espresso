#
# Copyright (C) 2013,2014,2015,2016 The ESPResSo project
#
# This file is part of ESPResSo.
#
# ESPResSo is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ESPResSo is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


from libcpp.map cimport map
from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp.memory cimport unique_ptr
from libcpp cimport bool

cdef extern from "script_interface/Parameter.hpp" namespace "ScriptInterface":
    cdef cppclass ParameterType:
        bool operator==(const ParameterType &a, const ParameterType &b)
        
cdef extern from "script_interface/Parameter.hpp" namespace "ScriptInterface::ParameterType":
    cdef ParameterType BOOL
    cdef ParameterType INT
    cdef ParameterType DOUBLE
    cdef ParameterType STRING
    cdef ParameterType INT_VECTOR
    cdef ParameterType DOUBLE_VECTOR

cdef extern from "script_interface/Parameter.hpp" namespace "ScriptInterface":
    cdef cppclass Parameter:
        ParameterType type()
        int n_elements()
        bool required()

cdef extern from "script_interface/ScriptInterface.hpp":
    cdef cppclass Variant:
        Variant()
        Variant(const Variant&)
        Variant &operator=(const Variant&)

cdef extern from "script_interface/ScriptInterface.hpp" namespace "boost":
    T get[T](const Variant &)

cdef extern from "script_interface/ScriptInterface.hpp" namespace "ScriptInterface":
    Variant make_variant[T](const T& x)
    
    cdef cppclass ScriptInterfaceBase:
        const string name()
        map[string, Variant] get_parameters()
        map[string, Parameter] all_parameters() 
        Variant get_parameter(const string &name)
        void set_parameter(const string &name, const Variant &value)
        void set_parameters(map[string, Variant] &parameters)
    
cdef extern from "utils/Factory.hpp" namespace "Utils":
    unique_ptr[T] factory_make[T](const string& name)

cdef class PScriptInterface:
    cdef ScriptInterfaceBase *sip
    cdef map[string, Parameter] parameters
    cdef Variant make_variant(self, ParameterType type, value)

