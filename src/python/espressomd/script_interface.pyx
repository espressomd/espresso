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
from .utils import to_char_pointer, to_str
from .utils cimport Vector3d, make_array_locked, handle_errors

cdef class PObjectId:
    """Python interface to a core ObjectId object."""

    cdef ObjectId id

    def __richcmp__(PObjectId a, PObjectId b, op):
        if op == 2:
            return a.id == b.id
        else:
            raise NotImplementedError

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
    name : :obj:`str`
        Name of the core class to instantiate (method 1).
    \*\*kwargs
        Parameters for the core class constructor (method 1).
    oid : :class:`PObjectId`
        Object id of an existing core object (method 2).
    policy : :obj:`str`, \{'GLOBAL', 'LOCAL'\}
        Creation policy.

    Attributes
    ----------

    sip: shared_ptr
        Pointer to a ScriptInterface object in the core.
    policy_: :obj:`str`
        Creation policy.

    """

    def __init__(self, name=None, policy="GLOBAL", oid=None, **kwargs):
        cdef CreationPolicy policy_

        if policy == "GLOBAL":
            policy_ = GLOBAL
        elif policy == "LOCAL":
            policy_ = LOCAL
        else:
            raise Exception("Unknown policy '{}'.".format(policy))

        if oid:
            self.set_sip_via_oid(oid)
        else:
            self.set_sip(make_shared(to_char_pointer(name), policy_))
            self.sip.get().construct(self._sanitize_params(kwargs))

    def __richcmp__(a, b, op):
        if op == 2:
            return a.id() == b.id()
        else:
            raise NotImplementedError

    def _ref_count(self):
        return self.sip.use_count()

    def _valid_parameters(self):
        return [to_str(p.data()) for p in self.sip.get().valid_parameters()]

    cdef set_sip(self, shared_ptr[ScriptInterfaceBase] sip):
        self.sip = sip

    def set_sip_via_oid(self, PObjectId id):
        """Set the shared_ptr to an existing core ScriptInterface object via
        its object id.
        """
        oid = id.id
        try:
            ptr = get_instance(oid).lock()
            self.set_sip(ptr)
        except BaseException:
            raise Exception("Could not get sip for given_id")

    def id(self):
        """Return the core class object id (:class:`PObjectId`)."""
        oid = PObjectId()
        oid.id = self.sip.get().id()
        return oid

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
        return to_str(self.sip.get().name())

    def _serialize(self):
        return self.sip.get().serialize()

    def _unserialize(self, state):
        cdef shared_ptr[ScriptInterfaceBase] so_ptr = ScriptInterfaceBase.unserialize(state)
        self.set_sip(so_ptr)

    cdef VariantMap _sanitize_params(self, in_params) except *:
        cdef VariantMap out_params
        cdef Variant v

        valid_params = self._valid_parameters()

        for pname in in_params:
            if pname in valid_params:
                out_params[to_char_pointer(pname)] = python_object_to_variant(
                    in_params[pname])
            else:
                raise RuntimeError("Unknown parameter '{}'".format(pname))

        return out_params

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

cdef Variant python_object_to_variant(value):
    """Convert Python objects to C++ Variant objects."""
    cdef Variant v
    cdef vector[Variant] vec
    cdef PObjectId oid

    if value is None:
        return Variant()

    # The order is important, the object character should be preserved
    # even if the PScriptInterface derived class is iterable.
    if isinstance(value, PScriptInterface):
        # Map python object to id
        oid = value.id()
        return make_variant[ObjectId](oid.id)
    elif hasattr(value, '__iter__') and not(type(value) == str):
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
        raise TypeError("Unknown type for conversion to Variant")

cdef variant_to_python_object(const Variant & value) except +:
    """Convert C++ Variant objects to Python objects."""
    cdef ObjectId oid
    cdef vector[Variant] vec
    cdef shared_ptr[ScriptInterfaceBase] ptr
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
    if is_type[Vector3d](value):
        return make_array_locked(get_value[Vector3d](value))
    if is_type[ObjectId](value):
        # Get the id and build a corresponding object
        oid = get_value[ObjectId](value)

        # ObjectId is nullable, and the default id corresponds to "null".
        if oid != ObjectId():
            ptr = get_instance(oid).lock()

            if not ptr:
                raise Exception("Object failed to exist.")

            so_name = to_str(ptr.get().name())
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

            poid = PObjectId()
            poid.id = ptr.get().id()

            return pclass(oid=poid)
        else:
            return None
    if is_type[vector[Variant]](value):
        vec = get_value[vector[Variant]](value)
        res = []

        for i in vec:
            res.append(variant_to_python_object(i))

        return res

    raise TypeError("Unknown type")


def _unpickle_so_class(so_name, state):
    cdef shared_ptr[ScriptInterfaceBase] sip = ScriptInterfaceBase.unserialize(state)

    poid = PObjectId()
    poid.id = sip.get().id()

    so = _python_class_by_so_name[so_name](oid=poid)
    so.define_bound_methods()

    return so


class ScriptInterfaceHelper(PScriptInterface):

    """
    Base class from which to derive most interfaces to core ScriptInterface
    classes.
    """

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
            "Class " + self.__class__.__name__ + " does not have an attribute " + attr)

    def __setattr__(self, attr, value):
        if attr in self._valid_parameters():
            self.set_params(**{attr: value})
        else:
            self.__dict__[attr] = value

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
    :class:`~espressomd.constraints.Constraints` and
    :class:`~espressomd.lbboundaries.LBBoundaries`. Derived classes must
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

    def __reduce__(self):
        res = []
        for item in self.__iter__():
            res.append(item)

        return (_unpickle_script_object_registry,
                (self._so_name, self.get_params(), res))


def _unpickle_script_object_registry(so_name, params, items):
    so = _python_class_by_so_name[so_name](**params)
    for item in items:
        so.add(item)
    return so


# Map from script object names to their corresponding python classes
_python_class_by_so_name = {}


def script_interface_register(c):
    """
    Decorator used to register script interface classes.
    This will store a name-to-class relationship in a registry, so that
    parameters of type object can be instantiated as the correct python class
    """
    if not hasattr(c, "_so_name"):
        raise Exception("Python classes representing a script object must "
                        "define an _so_name attribute at class level")
    _python_class_by_so_name[c._so_name] = c
    return c


def init():
    initialize()
