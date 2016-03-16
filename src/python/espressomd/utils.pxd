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
import numpy as np
cimport numpy as np

cdef extern from "stdlib.h":
    void free(void * ptr)
    void * malloc(size_t size)
    void * realloc(void * ptr, size_t size)

cdef extern from "utils.hpp":
    ctypedef struct int_list "IntList":
        int * e
        int n

    ctypedef struct double_list "DoubleList":
        double * e
        int n

    cdef void init_intlist(int_list * il)
    cdef void alloc_intlist(int_list * il, int size)
    cdef void realloc_intlist(int_list * il, int size)

cdef int_list * create_int_list_from_python_object(obj)
cdef np.ndarray create_nparray_from_int_list(int_list * il)
cdef np.ndarray create_nparray_from_double_list(double_list * dl)
cdef np.ndarray create_nparray_from_double_array(double * x, int n)
cdef check_type_or_throw_except(x, n, t, msg)
cdef check_range_or_except(x, v_min, incl_min, v_max, incl_max)
