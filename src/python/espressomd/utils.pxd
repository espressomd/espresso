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
#
from __future__ import print_function, absolute_import
import numpy as np
cimport numpy as np

from libcpp.string cimport string  # import std::string as string
from libcpp.vector cimport vector  # import std::vector as vector

cdef extern from "stdlib.h":
    void free(void * ptr)
    void * malloc(size_t size)
    void * realloc(void * ptr, size_t size)

cdef extern from "utils/List.hpp":
    cppclass int_list "IntList":
        int_list()
        int_list(int)
        int_list(int, int)

        int& operator[](int)
        void resize(int)
        void push_back(int)

        int * e
        unsigned n

    cppclass double_list "DoubleList":
        double_list()
        double_list(int)
        double_list(int, double)

        double& operator[](int)

        double * e
        unsigned n

cdef extern from "utils/Histogram.hpp" namespace "Utils":
    cdef void unravel_index(const int* const len_dims, const int ndims, const int flattened_index, int* unravelled_index_out)

cdef int_list create_int_list_from_python_object(obj)
cdef np.ndarray create_nparray_from_int_list(int_list * il)
cdef np.ndarray create_nparray_from_double_list(double_list * dl)
cdef np.ndarray create_nparray_from_double_array(double * x, int n)
cdef check_type_or_throw_except(x, n, t, msg)
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
