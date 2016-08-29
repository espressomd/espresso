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
from libcpp.memory cimport shared_ptr
from libcpp.memory cimport weak_ptr
from libcpp cimport bool

cdef extern from "Vector.hpp":
    cdef cppclass Vector3d:
        Vector3d()
        Vector3d(vector[double])
        vector[double] as_vector()
    cdef cppclass Vector2d:
        Vector2d()
        Vector2d(vector[double])
        vector[double] as_vector()

cdef extern from "script_interface/Parameter.hpp" namespace "ScriptInterface":
    cdef cppclass ParameterType:
        bool operator == (const ParameterType & a, const ParameterType & b)

cdef extern from "script_interface/Parameter.hpp" namespace "ScriptInterface::ParameterType":
    cdef ParameterType BOOL
    cdef ParameterType INT
    cdef ParameterType DOUBLE
    cdef ParameterType STRING
    cdef ParameterType INT_VECTOR
    cdef ParameterType DOUBLE_VECTOR
    cdef ParameterType VECTOR3D
    cdef ParameterType VECTOR2D
    cdef ParameterType OBJECT

cdef extern from "script_interface/Parameter.hpp" namespace "ScriptInterface":
    cdef cppclass Parameter:
        ParameterType type()
        int n_elements()
        bool required()

cdef extern from "script_interface/ScriptInterface.hpp" namespace "ScriptInterface":
    void initialize()
    cdef cppclass Variant:
        Variant()
        Variant(const Variant & )
        Variant & operator = (const Variant &)
        int which()

cdef extern from "script_interface/ScriptInterface.hpp" namespace "boost":
    T get[T](const Variant &) except +

cdef extern from "script_interface/ScriptInterface.hpp" namespace "ScriptInterface":
    cdef cppclass OId:
        OId(int)
        int id

    Variant make_variant[T](const T & x)

    cdef cppclass ScriptInterfaceBase:
        const string name()
        map[string, Variant] get_parameters()
        map[string, Parameter] all_parameters()
        Variant get_parameter(const string & name)
        void set_parameter(const string & name, const Variant & value)
        void set_parameters(map[string, Variant] & parameters)
        Variant call_method(const string & name, const map[string, Variant] & parameters)
        int id()

cdef extern from "script_interface/ScriptInterface.hpp" namespace "ScriptInterface::ScriptInterfaceBase":
    shared_ptr[ScriptInterfaceBase] make_shared(const string & name)
    weak_ptr[ScriptInterfaceBase] get_instance(int id)

cdef class PScriptInterface:
    cdef shared_ptr[ScriptInterfaceBase] sip
    cdef map[string, Parameter] parameters
    cdef Variant make_variant(self, ParameterType type, value)
    cdef set_sip(self, shared_ptr[ScriptInterfaceBase] sip)
    cdef variant_to_python_object(self, Variant value)
