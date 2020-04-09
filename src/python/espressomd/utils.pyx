#
# Copyright (C) 2013-2019 The ESPResSo project
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
cimport cython
import numpy as np
from libcpp.vector cimport vector
include "myconfig.pxi"

cpdef check_type_or_throw_except(x, n, t, msg):
    """
    Checks that x is of type t and that n values are given, otherwise throws
    ValueError with the message msg. If x is an array/list/tuple, the type
    checking is done on the elements, and all elements are checked. Integers
    are accepted when a float was asked for.

    """
    # Check whether x is an array/list/tuple or a single value
    if n > 1:
        if hasattr(x, "__getitem__"):
            for i in range(len(x)):
                if not isinstance(x[i], t):
                    if not ((t == float and is_valid_type(x[i], int))
                            or (t == float and issubclass(type(x[i]), np.integer))) \
                            and not (t == int and issubclass(type(x[i]), np.integer)):
                        raise ValueError(
                            msg + " -- Item " + str(i) + " was of type " + type(x[i]).__name__)
        else:
            # if n>1, but the user passed a single value, also throw exception
            raise ValueError(
                msg + " -- A single value was given but " + str(n) + " were expected.")
    else:
        # N=1 and a single value
        if not isinstance(x, t):
            if not (t == float and is_valid_type(x, int)) and not (
                    t == int and issubclass(type(x), np.integer)):
                raise ValueError(msg + " -- Got an " + type(x).__name__)


cdef np.ndarray create_nparray_from_double_array(double * x, int len_x):
    """
    Returns a numpy array from double array

    Parameters
    ----------
    x : double* which is to be converted
    len_x: len of array

    """
    numpyArray = np.zeros(len_x)
    for i in range(len_x):
        numpyArray[i] = x[i]
    return numpyArray

cdef check_range_or_except(D, name, v_min, incl_min, v_max, incl_max):
    """
    Checks that x is in range [v_min,v_max] (include boundaries via
    incl_min/incl_max = true) or throws a ValueError. v_min/v_max = 'inf' to
    disable limit.

    """
    x = D[name]

    # Array/list/tuple
    if hasattr(x, "__len__"):
        if (v_min != "inf" and ((incl_min and not all(v >= v_min for v in x))
                                or (not incl_min and not all(v > v_min for v in x)))) or \
           (v_max != "inf" and ((incl_max and not all(v <= v_max for v in x))
                                or (not incl_max and not all(v < v_max for v in x)))):
            raise ValueError("In " + name + ": Some values in " + str(x) + "are out of range " +
                             ("[" if incl_min else "]") + str(v_min) + "," + str(v_max) + ("]" if incl_max else "["))
    # Single Value
    else:
        if (v_min != "inf" and ((incl_min and x < v_min) or (not incl_min and x <= v_min)) or
                v_max != "inf" and ((incl_max and x > v_max) or (not incl_max and x >= v_max))):
            raise ValueError("In " + name + ": Value " + str(x) + " is out of range " + ("[" if incl_min else "]") +
                             str(v_min) + "," + str(v_max) + ("]" if incl_max else "["))


def to_char_pointer(s):
    """
    Returns a char pointer which contains the information of the provided python string.

    Parameters
    ----------
    s : :obj:`str`

    """
    if isinstance(s, unicode):
        s = ( < unicode > s).encode('utf8')
    return s


def to_str(s):
    """
    Returns a python string.

    Parameters
    ----------
    s : char*

    """
    if isinstance(s, unicode):
        return < unicode > s
    elif isinstance(s, bytes):
        return ( < bytes > s).decode('ascii')
    else:
        raise ValueError('Unknown string type {}'.format(type(s)))


class array_locked(np.ndarray):

    """
    Returns a non-writable numpy.ndarray with a special error message upon usage
    of __setitem__  or in-place operators. Cast return in __get__ of array
    properties to array_locked to prevent these operations.

    """

    ERR_MSG = "ESPResSo array properties return non-writable arrays \
and can only be modified as a whole, not in-place or component-wise. \
Use numpy.copy(<ESPResSo array property>) to get a writable copy."

    def __new__(cls, input_array):
        obj = np.asarray(input_array).view(cls)
        obj.flags.writeable = False
        return obj

    def __add__(self, other):
        return np.copy(self) + other

    def __radd__(self, other):
        return other + np.copy(self)

    def __sub__(self, other):
        return np.copy(self) - other

    def __rsub__(self, other):
        return other - np.copy(self)

    def __repr__(self):
        return repr(np.array(self))

    def __setitem__(self, i, val):
        raise ValueError(array_locked.ERR_MSG)

    def __iadd__(self, val):
        raise ValueError(array_locked.ERR_MSG)

    def __isub__(self, val):
        raise ValueError(array_locked.ERR_MSG)

    def __imul__(self, val):
        raise ValueError(array_locked.ERR_MSG)

    def __idiv__(self, val):
        raise ValueError(array_locked.ERR_MSG)

    def __itruediv__(self, val):
        raise ValueError(array_locked.ERR_MSG)

    def __ifloordiv__(self, val):
        raise ValueError(array_locked.ERR_MSG)

    def __imod__(self, val):
        raise ValueError(array_locked.ERR_MSG)

    def __ipow__(self, val):
        raise ValueError(array_locked.ERR_MSG)

    def __ilshift__(self, val):
        raise ValueError(array_locked.ERR_MSG)

    def __irshift__(self, val):
        raise ValueError(array_locked.ERR_MSG)

    def __iand__(self, val):
        raise ValueError(array_locked.ERR_MSG)

    def __ior__(self, val):
        raise ValueError(array_locked.ERR_MSG)

    def __ixor__(self, val):
        raise ValueError(array_locked.ERR_MSG)


cdef make_array_locked(Vector3d v):
    return array_locked([v[0], v[1], v[2]])


cdef Vector3d make_Vector3d(a):
    cdef Vector3d v
    for i, ai in enumerate(a):
        v[i] = ai
    return v


cpdef handle_errors(msg):
    """
    Gathers runtime errors.

    Parameters
    ----------
    msg: :obj:`str`
         Error message that is to be raised.

    """
    errors = mpi_gather_runtime_errors()
    # print all errors and warnings
    for err in errors:
        err.print()

    # raise an exception with the first error
    for err in errors:
        # Cast because cython does not support typed enums completely
        if < int > err.level() == < int > ERROR:
            raise Exception("{}: {}".format(msg, to_str(err.format())))


def nesting_level(obj):
    """
    Returns the maximal nesting level of an object.

    """

    if not isinstance(obj, (list, tuple, np.ndarray)):
        return 0

    obj = list(obj)

    max_level = 0
    for item in obj:
        max_level = max(max_level, nesting_level(item))

    return max_level + 1


def is_valid_type(value, t):
    """
    Extended checks for numpy int and float types.

    """
    if value is None:
        return False
    if t == int:
        return isinstance(value, (int, np.integer, np.long))
    elif t == float:
        if hasattr(np, 'float128'):
            return isinstance(
                value, (float, np.float16, np.float32, np.float64, np.float128, np.longdouble))
        return isinstance(
            value, (float, np.float16, np.float32, np.float64, np.longdouble))
    else:
        return isinstance(value, t)


def requires_experimental_features(reason):
    """Class decorator which makes instantiation conditional on
    ``EXPERIMENTAL_FEATURES`` being defined in myconfig.hpp."""

    def exception_raiser(self, *args, **kwargs):
        raise Exception(
            "Class " +
            self.__class__.__name__ +
            " is experimental. Define EXPERIMENTAL_FEATURES in myconfig.hpp to use it.\nReason: " +
            reason)

    def modifier(cls):
        cls.__init__ = exception_raiser
        return cls

    IF not EXPERIMENTAL_FEATURES:
        return modifier
    ELSE:
        # Return original class
        return lambda x: x
