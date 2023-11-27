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
#

from libcpp.string cimport string  # import std::string as string
from libcpp.vector cimport vector  # import std::vector as vector
from libcpp cimport bool as cbool

cpdef check_type_or_throw_except(x, n, t, msg)

cdef extern from "error_handling/RuntimeError.hpp" namespace "ErrorHandling::RuntimeError":
    cdef cppclass ErrorLevel:
        pass

cdef extern from "error_handling/RuntimeError.hpp" namespace "ErrorHandling::RuntimeError::ErrorLevel":
    cdef ErrorLevel WARNING
    cdef ErrorLevel ERROR

cdef extern from "error_handling/RuntimeError.hpp" namespace "ErrorHandling":
    cdef cppclass RuntimeError:
        string format()
        void print()
        ErrorLevel level()

cdef extern from "errorhandling.hpp" namespace "ErrorHandling":
    cdef vector[RuntimeError] mpi_gather_runtime_errors()

cpdef handle_errors(msg)

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

cdef make_array_locked(Vector3d)
cdef Vector3d make_Vector3d(a)
