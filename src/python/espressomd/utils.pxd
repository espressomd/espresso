#
# Copyright (C) 2013-2026 The ESPResSo project
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
#

from libcpp.string cimport string  # import std::string as string
from libcpp.vector cimport vector  # import std::vector as vector
from libcpp cimport bool as cbool

cdef extern from "<filesystem>" namespace "std::filesystem" nogil:
    cdef cppclass path:
        ctypedef char value_type
        path() except +
        path(const string & source) except +
        path & assign(const string & source) except +
        string native_string "string"() except +
        string generic_string() except +

cdef extern from "error_handling/RuntimeError.hpp" namespace "ErrorHandling::RuntimeError":
    cdef enum class ErrorLevel:
        WARNING
        ERROR

cdef extern from "error_handling/RuntimeError.hpp":
    cdef cppclass CoreRuntimeError "ErrorHandling::RuntimeError":
        string format()
        ErrorLevel level()

cdef extern from "errorhandling.hpp" namespace "ErrorHandling":
    cdef vector[CoreRuntimeError] mpi_gather_runtime_errors()
    cdef cbool mpi_poll_runtime_messages()

cdef extern from "utils/Vector.hpp" namespace "Utils":
    cppclass Vector2d:
        double & operator[](int i)
        double * data()

    cppclass Vector3d:
        double & operator[](int i)
        double * data()

    cppclass Vector4d:
        double & operator[](int i)
        double * data()

    cppclass Vector3b:
        cbool & operator[](int i)
        cbool * data()

    cppclass Vector3i:
        int & operator[](int i)
        int * data()
