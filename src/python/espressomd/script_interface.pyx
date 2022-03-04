# Copyright (C) 2010-2019 The ESPResSo project
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
import numpy as np
from .utils import to_char_pointer, to_str, handle_errors, array_locked, is_valid_type
from .utils cimport Vector2d, Vector3d, Vector4d, make_array_locked
cimport cpython.object

from libcpp.memory cimport make_shared

cdef shared_ptr[ContextManager] _om

cdef class PObjectRef:
    def __richcmp__(PObjectRef a, PObjectRef b, int op):
        if op == cpython.object.Py_EQ:
            return a.sip == b.sip
        elif op == cpython.object.Py_NE:
            return a.sip != b.sip
        else:
            raise NotImplementedError

    cdef shared_ptr[ObjectHandle] sip

    def print_sip(self):
        print( < long > (self.sip.get()))

cdef class PScriptInterface:

    """
    Python interface to a core ScriptInterface object. The core ScriptInterface
    class is itself an interface for other core classes, such as constraints,
    shapes, observables, etc.

    This class can be instantiated in two ways: (1) with the object id of an
    existing core ScriptInterface object, or (2) with parameters to construct
    a new ScriptInterface object in the core.

    Parameters
    ----------
    sip : :class:`PObjectRef`
        Object id of an existing core object (method 1).
    name : :obj:`str`
        Name of the core class to instantiate (method 2).
    \*\*kwargs
        Parameters for the core class constructor (method 2).
    policy : :obj:`str`, \{'GLOBAL', 'LOCAL'\}
        Creation policy. The managed object exists either on all MPI nodes
        with 'GLOBAL' (default), or only on the head node with 'LOCAL'.

    Attributes
    ----------

    sip: :class:`PObjectRef`
        Pointer to a ScriptInterface object in the core.

    """

    cdef shared_ptr[ObjectHandle] sip
    cdef set_sip(self, shared_ptr[ObjectHandle] sip)

    def __init__(self, name=None, policy="GLOBAL", sip=None, **kwargs):
        cdef CreationPolicy policy_
        cdef PObjectRef sip_
        cdef VariantMap out_params

        if policy == "GLOBAL":
            policy_ = GLOBAL
        elif policy == "LOCAL":
            policy_ = LOCAL
        else:
            raise Exception(f"Unknown policy '{policy}'.")

        if sip:
            sip_ = sip
            self.sip = sip_.sip
        else:
            global _om
            for pname in kwargs:
                out_params[to_char_pointer(pname)] = python_object_to_variant(
                    kwargs[pname])
            self.set_sip(
                _om.get().make_shared(
                    policy_,
                    to_char_pointer(name),
                    out_params))

    def __richcmp__(a, b, op):
        cls = PScriptInterface
        are_equality_comparable = isinstance(a, cls) and isinstance(b, cls)
        are_equal = are_equality_comparable and (a.get_sip() == b.get_sip())
        if op == cpython.object.Py_EQ:
            return are_equal
        elif op == cpython.object.Py_NE:
            return not are_equal
        else:
            raise NotImplementedError

    def _ref_count(self):
        return self.sip.use_count()

    def _valid_parameters(self):
        return [to_str(p.data()) for p in self.sip.get().valid_parameters()]

    def get_sip(self):
        """
        Get pointer to the core object.
        """

        ret = PObjectRef()
        ret.sip = self.sip

        return ret

    cdef set_sip(self, shared_ptr[ObjectHandle] sip):
        """
        Set the shared_ptr to an existing core object.
        """

        self.sip = sip

    def call_method(self, method, **kwargs):
        """
        Call a method of the core class.

        Parameters
        ----------
        method : Creation policy.
            Name of the core method.
        \*\*kwargs
            Arguments for the method.
        """
        cdef VariantMap parameters

        for name in kwargs:
            parameters[to_char_pointer(name)] = python_object_to_variant(
                kwargs[name])

        res = variant_to_python_object(
            self.sip.get().call_method(to_char_pointer(method), parameters))
        handle_errors("")
        return res

    def name(self):
        """Return name of the core class."""
        return to_str(self.sip.get().name().data())

    def _serialize(self):
        global _om
        return _om.get().serialize(self.sip.get())

    def _unserialize(self, state):
        global _om
        cdef shared_ptr[ObjectHandle] so_ptr = _om.get().deserialize(state)
        self.set_sip(so_ptr)

    def set_params(self, **kwargs):
        for name, value in kwargs.items():
            self.sip.get().set_parameter(to_char_pointer(name),
                                         python_object_to_variant(value))

    def get_parameter(self, name):
        cdef Variant value = self.sip.get().get_parameter(to_char_pointer(name))
        return variant_to_python_object(value)

    def get_params(self):
        odict = {}

        for pair in self.sip.get().get_parameters():
            odict[to_str(pair.first)] = variant_to_python_object(
                pair.second)

        return odict


class array_variant(np.ndarray):

    """
    Returns a numpy.ndarray that will be serialized as a ``std::vector``.

    """

    def __new__(cls, input_array):
        obj = np.asarray(input_array).view(cls)
        return obj


cdef Variant python_object_to_variant(value) except +:
    """Convert Python objects to C++ Variant objects."""

    cdef vector[Variant] vec
    cdef vector[int] vec_int
    cdef vector[double] vec_double
    cdef unordered_map[int, Variant] vmap
    cdef PObjectRef oref
    cdef int[::1] view_int
    cdef int * data_int
    cdef double[::1] view_double
    cdef double * data_double

    if value is None:
        return Variant()

    # The order is important, the object character should
    # be preserved even if the PScriptInterface derived class
    # is iterable.
    if isinstance(value, PScriptInterface):
        oref = value.get_sip()
        return make_variant(oref.sip)
    elif isinstance(value, dict):
        for k, v in value.items():
            if not isinstance(k, int):
                raise TypeError(
                    f"No conversion from type dict_item([({type(k).__name__}, {type(v).__name__})]) to Variant[std::unordered_map<int, Variant>]")
            vmap[k] = python_object_to_variant(v)
        return make_variant[unordered_map[int, Variant]](vmap)
    elif isinstance(value, array_variant) and np.issubdtype(value.dtype, np.signedinteger):
        view_int = np.ascontiguousarray(value, dtype=np.int32)
        data_int = &view_int[0]
        vec_int.assign(data_int, data_int + len(view_int))
        return make_variant[vector[int]](vec_int)
    elif isinstance(value, array_variant) and np.issubdtype(value.dtype, np.floating):
        view_double = np.ascontiguousarray(value, dtype=np.float64)
        data_double = &view_double[0]
        vec_double.assign(data_double, data_double + len(view_double))
        return make_variant[vector[double]](vec_double)
    elif hasattr(value, '__iter__') and type(value) != str:
        for e in value:
            vec.push_back(python_object_to_variant(e))
        return make_variant[vector[Variant]](vec)
    elif type(value) == str:
        return make_variant[string](to_char_pointer(value))
    elif type(value) == type(True):
        return make_variant[bool](value)
    elif np.issubdtype(np.dtype(type(value)), np.signedinteger):
        return make_variant[int](value)
    elif np.issubdtype(np.dtype(type(value)), np.floating):
        return make_variant[double](value)
    else:
        raise TypeError(
            f"No conversion from type {type(value).__name__} to Variant")

cdef variant_to_python_object(const Variant & value) except +:
    """Convert C++ Variant objects to Python objects."""

    cdef vector[Variant] vec
    cdef unordered_map[int, Variant] vmap
    cdef shared_ptr[ObjectHandle] ptr
    cdef Vector2d vec2d
    cdef Vector4d vec4d
    if is_none(value):
        return None
    if is_type[bool](value):
        return get_value[bool](value)
    if is_type[int](value):
        return get_value[int](value)
    if is_type[double](value):
        return get_value[double](value)
    if is_type[string](value):
        return to_str(get_value[string](value))
    if is_type[vector[int]](value):
        return get_value[vector[int]](value)
    if is_type[vector[double]](value):
        return get_value[vector[double]](value)
    if is_type[Vector4d](value):
        vec4d = get_value[Vector4d](value)
        return array_locked([vec4d[0], vec4d[1], vec4d[2], vec4d[3]])
    if is_type[Vector3d](value):
        return make_array_locked(get_value[Vector3d](value))
    if is_type[Vector2d](value):
        vec2d = get_value[Vector2d](value)
        return array_locked([vec2d[0], vec2d[1]])
    if is_type[shared_ptr[ObjectHandle]](value):
        # Get the id and build a corresponding object
        ptr = get_value[shared_ptr[ObjectHandle]](value)

        if ptr:
            so_name = to_str(ptr.get().name().data())
            if not so_name:
                raise Exception(
                    "Script object without name returned from the core")

            # Look up python type for object
            try:
                pclass = _python_class_by_so_name[so_name]
            except KeyError:
                # Fallback class, if nothing more specific is registered
                # for the script object name
                pclass = ScriptInterfaceHelper

            pptr = PObjectRef()
            pptr.sip = ptr

            return pclass(sip=pptr)
        else:
            return None
    if is_type[vector[Variant]](value):
        vec = get_value[vector[Variant]](value)
        res = []

        for i in vec:
            res.append(variant_to_python_object(i))

        return res
    if is_type[unordered_map[int, Variant]](value):
        vmap = get_value[unordered_map[int, Variant]](value)
        res = {}

        for kv in vmap:
            res[kv.first] = variant_to_python_object(kv.second)

        return res

    if is_type[size_t](value):
        return get_value[size_t](value)

    raise TypeError("Unknown type")


def _unpickle_so_class(so_name, state):
    cdef PObjectRef so_ptr
    so_ptr = PObjectRef()
    global _om
    so_ptr.sip = _om.get().deserialize(state)

    so = _python_class_by_so_name[so_name](sip=so_ptr)
    so.define_bound_methods()

    return so


class ScriptInterfaceHelper(PScriptInterface):
    _so_name = None
    _so_bind_methods = ()
    _so_creation_policy = "GLOBAL"

    def __init__(self, **kwargs):
        super().__init__(self._so_name, policy=self._so_creation_policy,
                         **kwargs)
        self.define_bound_methods()

    def __reduce__(self):
        return (_unpickle_so_class, (self._so_name, self._serialize()))

    def __dir__(self):
        return self.__dict__.keys() + self._valid_parameters()

    def __getattr__(self, attr):
        if attr in self._valid_parameters():
            return self.get_parameter(attr)

        if attr in self.__dict__:
            return self.__dict__[attr]

        raise AttributeError(
            f"Object '{self.__class__.__name__}' has no attribute '{attr}'")

    def __setattr__(self, attr, value):
        if attr in self._valid_parameters():
            self.set_params(**{attr: value})
        else:
            super().__setattr__(attr, value)

    def __delattr__(self, attr):
        if attr in self._valid_parameters():
            raise RuntimeError(f"Parameter '{attr}' is read-only")
        else:
            super().__delattr__(attr)

    def generate_caller(self, method_name):
        def template_method(**kwargs):
            res = self.call_method(method_name, **kwargs)
            return res

        return template_method

    def define_bound_methods(self):
        for method_name in self._so_bind_methods:
            setattr(self, method_name, self.generate_caller(method_name))


class ScriptObjectRegistry(ScriptInterfaceHelper):
    """
    Base class for container-like classes such as
    :class:`~espressomd.constraints.Constraints`. Derived classes must
    implement an ``add()`` method which adds a single item to the container.

    The core class should derive from ScriptObjectRegistry or provide
    ``"get_elements"`` and ``"size"`` as callable methods.
    """

    def __getitem__(self, key):
        return self.call_method("get_elements")[key]

    def __iter__(self):
        elements = self.call_method("get_elements")
        for e in elements:
            yield e

    def __len__(self):
        return self.call_method("size")


class ScriptObjectMap(ScriptObjectRegistry):

    """
    Represents a script interface ObjectMap (dict)-like object

    Item acces via [key].
    """

    _key_type = int

    def __init__(self, *args, **kwargs):
        if args:
            params, (_unpickle_so_class, (_so_name, bytestring)) = args
            assert _so_name == self._so_name
            self = _unpickle_so_class(_so_name, bytestring)
            self.__setstate__(params)
        else:
            super().__init__(**kwargs)

    def remove(self, key):
        """
        Removes the element with the given key
        """
        # Validate key type
        if not is_valid_type(key, self._key_type):
            raise ValueError(
                f"Key has to be of type {self._key_type.__name__}")

        self.call_method("erase", key=key)

    def clear(self):
        """
        Remove all elements.

        """
        self.call_method("clear")

    def _get(self, key):
        if not is_valid_type(key, self._key_type):
            raise ValueError(
                f"Key has to be of type {self._key_type.__name__}")

        return self.call_method("get", key=key)

    def __getitem__(self, key):
        return self._get(key)

    def __setitem__(self, key, value):
        self.call_method("insert", key=key, object=value)

    def __delitem__(self, key):
        self.remove(key)

    def keys(self):
        return self.call_method("keys")

    def __iter__(self):
        for k in self.keys(): yield k

    def items(self):
        for k in self.keys(): yield k, self[k]

    @classmethod
    def _restore_object(cls, so_callback, so_callback_args, state):
        so = so_callback(*so_callback_args)
        so.__setstate__(state)
        return so

    def __reduce__(self):
        so_callback, (so_name, so_bytestring) = super().__reduce__()
        return (ScriptObjectMap._restore_object,
                (so_callback, (so_name, so_bytestring), self.__getstate__()))

    def __getstate__(self):
        return dict(self.items())

    def __setstate__(self, params):
        for key, val in params.items():
            self[key] = val


# Map from script object names to their corresponding python classes
_python_class_by_so_name = {}


def script_interface_register(c):
    """
    Decorator used to register script interface classes.
    This will store a name<->class relationship in a registry, so that
    parameters of type object can be instantiated as the correct python class.
    """

    if not hasattr(c, "_so_name"):
        raise Exception("Python classes representing a script object must "
                        "define an _so_name attribute at class level")

    _python_class_by_so_name[c._so_name] = c
    return c


cdef void init(MpiCallbacks & cb):
    cdef Factory[ObjectHandle] f

    initialize(& f)

    global _om
    _om = make_shared[ContextManager](cb, f)
