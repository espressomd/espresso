#
# Copyright (C) 2013-2018 The ESPResSo project
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
from __future__ import print_function, absolute_import
import numpy as np
cimport numpy as np

from libcpp.string cimport string  # import std::string as string
from libcpp.vector cimport vector  # import std::vector as vector

cdef extern from "stdlib.h":
    void free(void * ptr)

cdef extern from "utils/List.hpp" namespace "Utils":
    cppclass List[T]:
        List()
        List(size_t)
        List(size_t, const T & )

        T & operator[](size_t)
        void resize(size_t)
        void push_back(size_t)

        T * data()
        size_t size()

        T * e
        size_t n

cdef extern from "utils/Span.hpp" namespace "Utils":
    cppclass Span[T]:
        Span()
        Span(T *, size_t)

        T & operator[](size_t)

        T * begin()
        T * end()

        T * data()
        size_t size()

    Span[const T] make_const_span[T](T *, size_t)

cdef List[int] create_int_list_from_python_object(obj)
cdef np.ndarray create_nparray_from_int_list(const List[int] & il)
cdef np.ndarray create_nparray_from_double_array(double * x, int n)
cpdef check_type_or_throw_except(x, n, t, msg)
cdef check_range_or_except(D, x, v_min, incl_min, v_max, incl_max)

cdef extern from "RuntimeError.hpp" namespace "ErrorHandling::RuntimeError":
    cdef cppclass ErrorLevel:
        pass

cdef extern from "RuntimeError.hpp" namespace "ErrorHandling::RuntimeError::ErrorLevel":
    cdef ErrorLevel WARNING
    cdef ErrorLevel ERROR

cdef extern from "RuntimeError.hpp" namespace "ErrorHandling":
    cdef cppclass RuntimeError:
        string format()
        void print()
        ErrorLevel level()

cdef extern from "errorhandling.hpp" namespace "ErrorHandling":
    cdef vector[RuntimeError] mpi_gather_runtime_errors()

cpdef handle_errors(msg)

# https://github.com/cython/cython/blob/master/Cython/Includes/libcpp/limits.pxd
cdef extern from "<limits>" namespace "std" nogil:
    cdef cppclass numeric_limits[T]:
        @staticmethod
        T epsilon()

        @staticmethod
        T max()

cdef extern from "utils/Vector.hpp" namespace "Utils":
    cppclass Vector2d:
        pass
    cppclass Vector4d:
        pass

    cppclass Vector3i:
        int & operator[](int i)
        int * data()

    cppclass Vector3d:
        double & operator[](int i)
        double * data()
        Vector3d operator * (double i)
        Vector3d operator / (double i)

    cppclass Vector6d:
        double & operator[](int i)
        double * data()
        Vector6d operator * (double i)
        Vector6d operator / (double i)

    cppclass Vector19d:
        double & operator[](int i)
        double * data()

cdef make_array_locked(const Vector3d & v)
cdef Vector3d make_Vector3d(a)
