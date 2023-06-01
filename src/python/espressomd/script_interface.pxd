#
# Copyright (C) 2013-2022 The ESPResSo project
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


from libcpp.unordered_map cimport unordered_map
from libcpp.string cimport string
from libcpp.memory cimport shared_ptr
from libcpp.vector cimport vector
from libcpp cimport bool

from boost cimport string_ref

from .communication cimport MpiCallbacks

cdef extern from "utils/Factory.hpp" namespace "Utils":
    cdef cppclass Factory[T]:
        pass

cdef extern from "script_interface/ScriptInterface.hpp" namespace "ScriptInterface":
    cdef cppclass Variant:
        Variant()
        Variant(const Variant & )
        Variant & operator = (const Variant &)

    bool is_type[T](const Variant &)
    bool is_none(const Variant &)
    ctypedef unordered_map[string, Variant] VariantMap

    Variant make_variant[T](const T & x)

    cdef cppclass ObjectHandle:
        VariantMap get_parameters() except +
        vector[string] get_valid_parameters() except +
        Variant get_parameter(const string & name) except +
        void set_parameter(const string & name, const Variant & value) except +
        Variant call_method(const string & name, const VariantMap & parameters) except +
        Variant call_method_nogil "call_method"(const string & name, const VariantMap & parameters) nogil except +
        string_ref name()

cdef extern from "script_interface/ContextManager.hpp" namespace "ScriptInterface::ContextManager":
    cdef cppclass CreationPolicy:
        pass

cdef extern from "script_interface/ContextManager.hpp" namespace "ScriptInterface::ContextManager::CreationPolicy":
    CreationPolicy LOCAL
    CreationPolicy GLOBAL

cdef extern from "script_interface/ContextManager.hpp" namespace "ScriptInterface":
    cdef cppclass ContextManager:
        ContextManager(MpiCallbacks & , const Factory[ObjectHandle] & )
        shared_ptr[ObjectHandle] make_shared(CreationPolicy, const string &, const VariantMap) except +
        shared_ptr[ObjectHandle] deserialize(const string &) except +
        string serialize(const ObjectHandle *) except +

cdef extern from "script_interface/initialize.hpp" namespace "ScriptInterface":
    void initialize(Factory[ObjectHandle] *)

cdef extern from "script_interface/get_value.hpp" namespace "ScriptInterface":
    T get_value[T](const Variant T)

cdef void init(MpiCallbacks &)
