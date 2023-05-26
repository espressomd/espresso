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
cimport numpy as np
import numpy as np


cdef _check_type_or_throw_except_assertion(x, t):
    return isinstance(x, t) or (t == int and is_valid_type(x, int)) or (
        t == float and (is_valid_type(x, int) or is_valid_type(x, float))) or (
        t == bool and is_valid_type(x, bool))


cpdef check_array_type_or_throw_except(x, n, t, msg):
    """
    Check that ``x`` is of type ``t`` and that ``n`` values are given,
    otherwise raise a ``ValueError`` with message ``msg``.
    Integers are accepted when a float was asked for.

    """
    if not hasattr(x, "__getitem__"):
        raise ValueError(
            msg + f" -- A single value was given but {n} were expected.")
    if len(x) != n:
        raise ValueError(
            msg + f" -- {len(x)} values were given but {n} were expected.")
    if isinstance(x, np.ndarray):
        value = x.dtype.type()  # default-constructed value of the same type
        if not _check_type_or_throw_except_assertion(value, t):
            raise ValueError(
                msg + f" -- Array was of type {type(value).__name__}")
        return
    for i in range(len(x)):
        if not _check_type_or_throw_except_assertion(x[i], t):
            raise ValueError(
                msg + f" -- Item {i} was of type {type(x[i]).__name__}")


cpdef check_type_or_throw_except(x, n, t, msg):
    """
    Check that ``x`` is of type ``t`` and that ``n`` values are given,
    otherwise raise a ``ValueError`` with message ``msg``. If ``x`` is an
    array/list/tuple, the type checking is done on the elements, and all
    elements are checked. If ``n`` is 1, ``x`` is assumed to be a scalar.
    Integers are accepted when a float was asked for.

    """
    # Check whether x is an array/list/tuple or a single value
    if n > 1:
        check_array_type_or_throw_except(x, n, t, msg)
    else:
        if not _check_type_or_throw_except_assertion(x, t):
            raise ValueError(msg + f" -- Got an {type(x).__name__}")


def to_char_pointer(s):
    """
    Returns a Cython bytes object which contains the information of the provided
    Python string. Cython bytes objects implicitly cast to raw char pointers.

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
        raise ValueError(f'Unknown string type {type(s)}')


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
            raise Exception(f"{msg}: {to_str(err.format())}")


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
    Extended checks for numpy int, float and bool types.
    Handles 0-dimensional arrays.

    """
    if value is None:
        return False
    if isinstance(value, np.ndarray) and value.shape == ():
        value = value[()]
    if t == int:
        return isinstance(value, (int, np.integer))
    elif t == float:
        float_types = [
            float, np.float16, np.float32, np.float64, np.longdouble]
        if hasattr(np, 'float128'):
            float_types.append(np.float128)
        return isinstance(value, tuple(float_types))
    elif t == bool:
        return isinstance(value, (bool, np.bool_))
    else:
        return isinstance(value, t)


def check_required_keys(required_keys, obtained_keys):
    a = required_keys
    b = obtained_keys
    if not set(a).issubset(b):
        raise ValueError(
            "The following keys have to be given as keyword arguments: "
            f"{sorted(a)}, got {sorted(b)} (missing {sorted(a - b)})")


def check_valid_keys(valid_keys, obtained_keys):
    a = valid_keys
    b = obtained_keys
    if not set(b).issubset(a):
        raise ValueError(
            "Only the following keys can be given as keyword arguments: "
            f"{sorted(a)}, got {sorted(b)} (unknown {sorted(b - a)})")
