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
cimport numpy as np
import numpy as np

cdef extern from "stdlib.h":
    void free(void * ptr)
    void * malloc(size_t size)
    void * realloc(void * ptr, size_t size)

cdef np.ndarray create_nparray_from_int_list(int_list * il):
    numpyArray = np.zeros(il.n)
    for i in range(il.n):
        numpyArray[i] = il.e[i]

    return numpyArray

cdef np.ndarray create_nparray_from_double_list(double_list * dl):
    numpyArray = np.zeros(dl.n)
    for i in range(dl.n):
        numpyArray[i] = dl.e[i]

cdef int_list * create_int_list_from_python_object(obj):
    cdef int_list * il
    il = <int_list * > malloc(sizeof(int_list))
    init_intlist(il)

    alloc_intlist(il, len(obj))
    for i in range(len(obj)):
        il.e[i] = obj[i]
        print il.e[i]
    return il

cdef check_type_or_throw_except(x, n, t, msg):
    """Checks that x is of type t and that n values are given, otherwise throws ValueError with the message msg.
       If x is an array/list/tuple, the type checking is done on the elements, and
       all elements are checked.
       Integers are accepted when a float was asked for.
       """
    # Check whether x is an array/list/tuple or a single value
    if n > 1:
        if hasattr(x, "__getitem__"):
            for i in range(len(x)):
                if not isinstance(x[i], t):
                    if not (t == float and isinstance(x[i], int)):
                        raise ValueError(
                            msg + " -- Item " + str(i) + " was of type " + type(x[i]).__name__)
        else:
            # if n>1, but the user passed a single value, also throw exception
            raise ValueError(
                msg + " -- A single value was given but " + str(n) + " were expected.")
    else:
        # N=1 and a single value
        if not isinstance(x, t):
            if not (t == float and isinstance(x, int)):
                raise ValueError(msg + " -- Got an " + type(x).__name__)


cdef np.ndarray create_nparray_from_double_array(double * x, int n):
    numpyArray = np.zeros(n)
    for i in range(n):
        numpyArray[i] = x[i]
    return numpyArray

cdef check_range_or_except(x, v_min, incl_min, v_max, incl_max):
    """Checks that x is in range [v_min,v_max] (inlude boundaries via inlc_min/incl_max = true) or throws a ValueError. v_min/v_max = 'inf' to disable limit """

    # Array/list/tuple
    if hasattr(x, "__len__"):
        if (v_min != "inf" and ((incl_min and not all(v >= v_min for v in x))
                                or not all(v > v_min for v in x))) or \
           (v_max != "inf" and ((incl_max and not all(v <= v_max for v in x))
                                or not all(v < v_max for v in x))):
            raise ValueError("Some values in " + str(x) + "are out of range " +
                             "[" if incl_min else "]" + str(v_min) + "," + str(v_max) + "]" if incl_max else "[")
    # Single Value
    else:
        if (v_min != "inf" and ((incl_min and not x >= v_min)
                                or not x > v_min)) or \
           (v_max != "inf" and ((incl_max and not x <= v_max)
                                or not x < v_max)):
            raise ValueError("Value " + str(x) + "is out of range " + "[" if incl_min else "]" + str(
                v_min) + "," + str(v_max) + "]" if incl_max else "[")
