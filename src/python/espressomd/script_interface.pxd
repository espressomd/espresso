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

cdef extern from "script_interface/Parameter.hpp" namespace "ScriptInterface":
    cdef cppclass ParameterType:
        bool operator == (const ParameterType & a, const ParameterType & b)

cdef extern from "script_interface/Parameter.hpp" namespace "ScriptInterface::ParameterType":
    cdef ParameterType NONE
    cdef ParameterType BOOL
    cdef ParameterType INT
    cdef ParameterType DOUBLE
    cdef ParameterType STRING
    cdef ParameterType INT_VECTOR
    cdef ParameterType DOUBLE_VECTOR
    cdef ParameterType OBJECTID
    cdef ParameterType VECTOR

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
    void transform_vectors(Variant &)
    string get_type_label(const Variant &)
    string get_type_label(ParameterType)

cdef extern from "script_interface/ScriptInterface.hpp" namespace "boost":
    T get[T](const Variant &) except +

cdef extern from "script_interface/ScriptInterface.hpp" namespace "ScriptInterface":
    cdef cppclass ObjectId:
        string to_string()
        bool operator==(const ObjectId& rhs) 

    Variant make_variant[T](const T & x)

    cdef cppclass ScriptInterfaceBase:
        const string name()
        void construct(map[string, Variant]) except +
        map[string, Variant] get_parameters() except +
        map[string, Parameter] valid_parameters() except +
        Variant get_parameter(const string & name) except +
        void set_parameter(const string & name, const Variant & value) except +
        void set_parameters(map[string, Variant] & parameters) except +
        Variant call_method(const string & name, const map[string, Variant] & parameters) except +
        ObjectId id() except +
        void set_state(map[string, Variant]) except +
        map[string, Variant] get_state() except +
        string serialize() except +
        @staticmethod
        shared_ptr[ScriptInterfaceBase] unserialize(const string &state) except +

cdef extern from "script_interface/ScriptInterface.hpp" namespace "ScriptInterface::ScriptInterfaceBase":
    cdef cppclass CreationPolicy:
        pass
    shared_ptr[ScriptInterfaceBase] make_shared(const string & name, CreationPolicy policy) except +
    weak_ptr[ScriptInterfaceBase] get_instance(ObjectId id) except +

cdef extern from "script_interface/ScriptInterface.hpp" namespace "ScriptInterface::ScriptInterfaceBase::CreationPolicy":
    CreationPolicy LOCAL
    CreationPolicy GLOBAL

cdef class PScriptInterface:
    cdef shared_ptr[ScriptInterfaceBase] sip
    cdef map[string, Parameter] parameters
    cdef set_sip(self, shared_ptr[ScriptInterfaceBase] sip)
    cdef variant_to_python_object(self, Variant value)  except +
    cdef Variant python_object_to_variant(self, value)
    cdef map[string, Variant] _sanitize_params(self, in_params)

