# Copyright (C) 2010-2022 The ESPResSo project
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
from libcpp cimport bool as cbool
from libcpp.unordered_map cimport unordered_map
from libcpp.utility cimport pair
from libcpp.vector cimport vector
from libcpp.memory cimport shared_ptr, make_shared
import numpy as np
cimport numpy as cnp
import pathlib
from . import utils
from .utils cimport Vector3b, Vector3i, Vector2d, Vector3d, Vector4d
from .utils cimport path
cimport cpython.object
cnp.import_array()


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
    **kwargs
        Parameters for the core class constructor (method 2).
    policy : :obj:`str`, {'GLOBAL', 'LOCAL'}
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
                out_params[utils.to_bytes(pname)] = python_object_to_variant(
                    kwargs[pname])
            self.set_sip(
                _om.get().make_shared(
                    policy_,
                    utils.to_bytes(name),
                    out_params))
            utils.handle_errors(f"Raised during instantiation of '{name}'")

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
        cdef ObjectHandle * handle = self.sip.get()
        return [utils.to_str(p.data()) for p in handle.valid_parameters()]

    def _has_parameter(self, name):
        return self.sip.get().has_parameter(utils.to_bytes(name))

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

    def call_method(self, method, handle_errors_message=None,
                    with_nogil=False, **kwargs):
        """
        Call a method of the core class.

        Parameters
        ----------
        method : :obj:`str`
            Name of the core method.
        handle_errors_message : :obj:`str`, optional
            Custom error message for runtime errors raised in a MPI context.
        with_nogil : :obj:`bool`, optional
            Run the core method without the GIL if ``True``.
        **kwargs
            Arguments for the method.
        """
        cdef ObjectHandle * handle = self.sip.get()
        cdef VariantMap parameters
        cdef Variant result

        for name, value in kwargs.items():
            parameters[utils.to_bytes(name)] = python_object_to_variant(value)

        # the internal buffer of a cython bytestring object can be accessed as
        # a raw char pointer, but then the bytestring object must be kept alive
        method_name_bytes_counted_reference = utils.to_bytes(method)
        cdef char * method_name_char = method_name_bytes_counted_reference

        if with_nogil:
            with nogil:
                result = handle.call_method_nogil(method_name_char, parameters)
        else:
            result = handle.call_method(method_name_char, parameters)
        result_py = variant_to_python_object(result)
        if handle_errors_message is None:
            handle_errors_message = f"Raised while calling method {method}()"
        utils.handle_errors(handle_errors_message)
        return result_py

    def name(self):
        """Return name of the core class."""
        return utils.to_str(self.sip.get().name().data())

    def _serialize(self):
        global _om
        return _om.get().serialize(self.sip.get())

    def _unserialize(self, state):
        global _om
        cdef shared_ptr[ObjectHandle] so_ptr = _om.get().deserialize(state)
        self.set_sip(so_ptr)

    def set_params(self, **kwargs):
        for name, value in kwargs.items():
            self.sip.get().set_parameter(utils.to_bytes(name),
                                         python_object_to_variant(value))
            utils.handle_errors(f"while setting parameter '{name}'")

    def get_parameter(self, name):
        cdef Variant value = self.sip.get().get_parameter(utils.to_bytes(name))
        return variant_to_python_object(value)

    def get_params(self):
        odict = {}

        for pair in self.sip.get().get_parameters():
            odict[utils.to_str(pair.first)] = variant_to_python_object(
                pair.second)

        return odict


class array_variant(np.ndarray):

    """
    Returns a numpy.ndarray that will be serialized as a ``std::vector``.

    """

    def __new__(cls, input_array):
        obj = np.asarray(input_array).view(cls)
        return obj


def fast_tiling(value, n):
    """
    Repeat a value (scalar or vector) multiple times.
    Based on the value type, either a NumPy n-dimensional array is returned,
    or a list of references to the original object (if mutable),
    or a list of copies of the original object (if immutable).
    """
    if isinstance(value, np.ndarray) and value.ndim == 1 and len(value) != 0 and \
            any(np.issubdtype(value.dtype, x) for x in (np.floating, np.signedinteger)) or \
            isinstance(value, (list, tuple)) and len(value) != 0 and \
            all(map(lambda x: isinstance(x, (int, np.signedinteger, float, np.floating))
                    and not isinstance(x, type(True)), value)):
        return np.tile(value, (n, 1))
    if isinstance(value, (int, np.signedinteger, float, np.floating)) and \
            not isinstance(value, type(True)):
        return array_variant(np.repeat(value, n))
    return n * [value]


cdef Variant python_object_to_variant(value) except *:
    """Convert Python objects to C++ Variant objects."""

    cdef vector[Variant] vec_variant
    cdef vector[int] vec_int
    cdef vector[double] vec_double
    cdef unordered_map[int, Variant] map_int2var
    cdef unordered_map[string, Variant] map_str2var
    cdef PObjectRef oref
    cdef int[::1] view_int
    cdef int[:, ::1] view_int_2d
    cdef int * data_int
    cdef double[::1] view_double
    cdef double[:, ::1] view_double_2d
    cdef double * data_double
    cdef path fs_path
    cdef size_t index
    cdef size_t nrows
    cdef size_t bufsize
    cdef Vector3d vector3d
    cdef Vector3i vector3i

    if value is None:
        return Variant()

    if isinstance(value, np.ndarray) and value.ndim == 0:
        value = value.item()

    # The order is important, the object character should
    # be preserved even if the PScriptInterface derived class
    # is iterable.
    if isinstance(value, PScriptInterface):
        oref = value.get_sip()
        return make_variant(oref.sip)
    if isinstance(value, dict):
        if all(map(lambda x: isinstance(x, (int, np.integer)), value.keys())):
            for key, value in value.items():
                map_int2var[int(key)] = python_object_to_variant(value)
            return make_variant[unordered_map[int, Variant]](map_int2var)
        if all(map(lambda x: isinstance(x, (str, bytes)), value.keys())):
            for key, value in value.items():
                key_bytes = utils.to_bytes(key)
                map_str2var[key_bytes] = python_object_to_variant(value)
            return make_variant[unordered_map[string, Variant]](map_str2var)
        for k, v in value.items():
            if not isinstance(k, (str, bytes, int, np.integer)):
                raise TypeError(
                    f"No conversion from type "
                    f"'dict_item([({type(k).__name__}, {type(v).__name__})])'"
                    f" to 'Variant[std::unordered_map<int, Variant>]' or"
                    f" to 'Variant[std::unordered_map<std::string, Variant>]'")
        assert False, "dev note: a type is missing in the for loop above"
    if isinstance(value, (str, bytes)):
        return make_variant[string](utils.to_bytes(value))
    if isinstance(value, pathlib.Path):
        fs_path.assign(utils.to_bytes(str(value)))
        return make_variant[path](fs_path)
    if isinstance(value, np.ndarray):
        if isinstance(value, array_variant):
            if np.issubdtype(value.dtype, np.signedinteger):
                view_int = np.ascontiguousarray(value, dtype=np.int32)
                data_int = &view_int[0]
                vec_int.assign(data_int, data_int + len(view_int))
                return make_variant[vector[int]](vec_int)
            if np.issubdtype(value.dtype, np.floating):
                view_double = np.ascontiguousarray(value, dtype=np.float64)
                data_double = &view_double[0]
                vec_double.assign(data_double, data_double + len(view_double))
                return make_variant[vector[double]](vec_double)
        if value.ndim == 1:
            if np.issubdtype(value.dtype, np.floating):
                vec_double.reserve(len(value))
                for e in value:
                    vec_double.push_back(e)
                return make_variant[vector[double]](vec_double)
            if np.issubdtype(value.dtype, np.signedinteger):
                vec_int.reserve(len(value))
                for e in value:
                    vec_int.push_back(e)
                return make_variant[vector[int]](vec_int)
        if value.ndim == 2:
            if np.issubdtype(value.dtype, np.signedinteger):
                nrows = value.shape[0]
                bufsize = value.shape[1]
                vec_variant.reserve(nrows)
                vec_int.reserve(bufsize)
                view_int_2d = np.ascontiguousarray(value, dtype=np.int32)
                for index in range(nrows):
                    data_int = &view_int_2d[index, 0]
                    vec_int.assign(data_int, data_int + bufsize)
                    vec_variant.emplace_back(
                        make_variant[vector[int]](vec_int))
                return make_variant[vector[Variant]](vec_variant)
            if np.issubdtype(value.dtype, np.floating):
                nrows = value.shape[0]
                bufsize = value.shape[1]
                vec_variant.reserve(nrows)
                vec_double.reserve(bufsize)
                view_double_2d = np.ascontiguousarray(value, dtype=np.float64)
                for index in range(nrows):
                    data_double = &view_double_2d[index, 0]
                    vec_double.assign(data_double, data_double + bufsize)
                    vec_variant.emplace_back(
                        make_variant[vector[double]](vec_double))
                return make_variant[vector[Variant]](vec_variant)
    if hasattr(value, "__iter__"):
        bufsize = len(value)
        if bufsize == 0:
            return make_variant[vector[Variant]](vec_variant)
        if all(map(lambda x: isinstance(x, (float, np.floating)), value)):
            if bufsize == 3:
                for index in range(bufsize):
                    vector3d[index] = value[index]
                return make_variant[Vector3d](vector3d)
            vec_double.reserve(bufsize)
            for e in value:
                vec_double.push_back(e)
            return make_variant[vector[double]](vec_double)
        if all(map(lambda x: isinstance(x, (int, np.integer))
                   and not isinstance(x, type(True)), value)):
            if bufsize == 3:
                for index in range(bufsize):
                    vector3i[index] = value[index]
                return make_variant[Vector3i](vector3i)
            vec_int.reserve(bufsize)
            for e in value:
                vec_int.push_back(e)
            return make_variant[vector[int]](vec_int)
        vec_variant.reserve(bufsize)
        for e in value:
            vec_variant.emplace_back(python_object_to_variant(e))
        return make_variant[vector[Variant]](vec_variant)
    if isinstance(value, (type(True), np.bool_)):
        return make_variant[cbool](value)
    if np.issubdtype(np.dtype(type(value)), np.signedinteger):
        return make_variant[int](value)
    if np.issubdtype(np.dtype(type(value)), np.floating):
        return make_variant[double](value)
    raise TypeError(
        f"No conversion from type '{type(value).__name__}' to 'Variant'")


cdef variant_to_python_object(const Variant & value):
    """Convert C++ Variant objects to Python objects."""

    cdef vector[Variant] vec
    cdef unordered_map[int, Variant] map_int2var
    cdef unordered_map[string, Variant] map_str2var
    cdef pair[int, Variant] pair_int2var
    cdef pair[string, Variant] pair_str2var
    cdef shared_ptr[ObjectHandle] ptr
    cdef Vector3b vec3b
    cdef Vector3i vec3i
    cdef Vector2d vec2d
    cdef Vector3d vec3d
    cdef Vector4d vec4d
    cdef cnp.ndarray[cnp.float64_t, ndim = 2] arrayNvec3d
    cdef size_t index
    cdef size_t nrows
    if is_none(value):
        return None
    if is_type[cbool](value):
        return get_value[cbool](value)
    if is_type[int](value):
        return get_value[int](value)
    if is_type[double](value):
        return get_value[double](value)
    if is_type[string](value):
        return utils.to_str(get_value[string](value))
    if is_type[path](value):
        filepath = utils.to_str(get_value[path](value).generic_string())
        return pathlib.Path(filepath)
    if is_type[vector[int]](value):
        return get_value[vector[int]](value)
    if is_type[vector[double]](value):
        return np.array(get_value[vector[double]](value))
    if is_type[Vector3b](value):
        vec3b = get_value[Vector3b](value)
        return utils.array_locked([vec3b[0], vec3b[1], vec3b[2]])
    if is_type[Vector3i](value):
        vec3i = get_value[Vector3i](value)
        return utils.array_locked([vec3i[0], vec3i[1], vec3i[2]])
    if is_type[Vector4d](value):
        vec4d = get_value[Vector4d](value)
        return utils.array_locked([vec4d[0], vec4d[1], vec4d[2], vec4d[3]])
    if is_type[Vector3d](value):
        vec3d = get_value[Vector3d](value)
        return utils.array_locked([vec3d[0], vec3d[1], vec3d[2]])
    if is_type[Vector2d](value):
        vec2d = get_value[Vector2d](value)
        return utils.array_locked([vec2d[0], vec2d[1]])
    if is_type[shared_ptr[ObjectHandle]](value):
        # Get the id and build a corresponding object
        ptr = get_value[shared_ptr[ObjectHandle]](value)

        if ptr:
            so_name = utils.to_str(ptr.get().name().data())
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
        nrows = vec.size()

        if (nrows > 0) and is_type[Vector3d](vec[0]):
            arrayNvec3d = np.empty((nrows, 3), dtype=np.float64)
            mixed_types = False
            for index in range(nrows):
                if not is_type[Vector3d](vec[index]):
                    mixed_types = True
                    break
                vec3d = get_value[Vector3d](vec[index])
                arrayNvec3d[index, 0] = vec3d[0]
                arrayNvec3d[index, 1] = vec3d[1]
                arrayNvec3d[index, 2] = vec3d[2]
            if not mixed_types:
                return arrayNvec3d

        res = []
        for index in range(nrows):
            res.append(variant_to_python_object(vec[index]))
        return res

    if is_type[unordered_map[int, Variant]](value):
        map_int2var = get_value[unordered_map[int, Variant]](value)
        res = {}

        for pair_int2var in map_int2var:
            res[pair_int2var.first] = variant_to_python_object(
                pair_int2var.second)

        return res

    if is_type[unordered_map[string, Variant]](value):
        map_str2var = get_value[unordered_map[string, Variant]](value)
        res = {}

        for pair_str2var in map_str2var:
            res[utils.to_str(pair_str2var.first)] = variant_to_python_object(
                pair_str2var.second)

        return res

    if is_type[size_t](value):
        return get_value[size_t](value)

    raise TypeError("Unknown type")


def _unpickle_so_class(so_name, state):
    cdef PObjectRef so_ptr
    so_ptr = PObjectRef()
    global _om
    so_ptr.sip = _om.get().deserialize(state)

    assert so_name in _python_class_by_so_name, \
        f"C++ class '{so_name}' is not associated to any Python class " \
        "(hint: the corresponding 'import espressomd.*' may be missing)"
    so = _python_class_by_so_name[so_name](sip=so_ptr)
    so.define_bound_methods()

    return so


class ScriptInterfaceHelper(PScriptInterface):
    _so_name = None
    _so_features = ()
    _so_bind_methods = ()
    _so_checkpointable = True
    _so_creation_policy = "GLOBAL"

    def __init__(self, **kwargs):
        cdef vector[string] features_vec
        if self._so_features:
            for feature in self._so_features:
                features_vec.push_back(utils.to_bytes(feature))
            check_features(features_vec)
        super().__init__(self._so_name, policy=self._so_creation_policy,
                         **kwargs)
        self.define_bound_methods()

    def __reduce__(self):
        assert self._so_checkpointable, f"Class '{self.__class__.__name__}' doesn't support checkpointing"  # nopep8
        return (_unpickle_so_class, (self._so_name, self._serialize()))

    def __dir__(self):
        return list(self.__dict__.keys()) + self._valid_parameters()

    def __getattr__(self, attr):
        if self._has_parameter(attr):
            return self.get_parameter(attr)

        if attr in self.__dict__:
            return self.__dict__[attr]

        raise AttributeError(
            f"Object '{self.__class__.__name__}' has no attribute '{attr}'")

    def __setattr__(self, attr, value):
        if self._has_parameter(attr):
            self.set_params(**{attr: value})
        else:
            super().__setattr__(attr, value)

    def __delattr__(self, attr):
        if self._has_parameter(attr):
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


class ScriptObjectList(ScriptInterfaceHelper):
    """
    Base class for container-like classes such as
    :class:`~espressomd.constraints.Constraints`. Derived classes must
    implement an ``add()`` method which adds a single item to the container.

    The core objects must be managed by a container derived from
    ``ScriptInterface::ObjectList``.

    """

    def __getitem__(self, key):
        return self.call_method("get_elements")[key]

    def __iter__(self):
        elements = self.call_method("get_elements")
        for e in elements:
            yield e

    def __len__(self):
        return self.call_method("size")


class ScriptObjectMap(ScriptInterfaceHelper):
    """
    Base class for container-like classes such as
    :class:`~espressomd.interactions.BondedInteractions`. Derived classes must
    implement an ``add()`` method which adds a single item to the container.

    The core objects must be managed by a container derived from
    ``ScriptInterface::ObjectMap``.

    """

    def remove(self, key):
        """
        Remove the element with the given key.
        This is a no-op if the key does not exist.
        """
        self.__delitem__(key)

    def clear(self):
        """
        Remove all elements.

        """
        self.call_method("clear")

    def __len__(self):
        return self.call_method("size")

    def __getitem__(self, key):
        return self.call_method("get", key=key)

    def __setitem__(self, key, value):
        self.call_method("insert", key=key, object=value)

    def __delitem__(self, key):
        self.call_method("erase", key=key)

    def keys(self):
        return self.call_method("keys")

    def __iter__(self):
        for k in self.keys():
            yield k

    def items(self):
        for k in self.keys():
            yield k, self[k]


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


cdef void init(const shared_ptr[MpiCallbacks] & cb):
    cdef Factory[ObjectHandle] f

    initialize(& f)

    global _om
    _om = make_shared[ContextManager](cb, f)


cdef void deinit():
    global _om
    _om.reset()
