from espressomd.utils import to_char_pointer, to_str, handle_errors
import numpy as np
from espressomd.utils import is_valid_type, array_locked
from espressomd.utils cimport Vector3d, make_array_locked

cdef class PObjectId(object):
    cdef ObjectId id

    def __richcmp__(PObjectId a, PObjectId b, op):
        if op == 2:
            return a.id == b.id
        else:
            raise NotImplementedError

cdef class PScriptInterface(object):
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
        """Set the shared_ptr to the script object in the core via the object id"""
        oid = id.id
        try:
            ptr = get_instance(oid).lock()
            self.set_sip(ptr)
        except:
            raise Exception("Could not get sip for given_id")

    def id(self):
        oid = PObjectId()
        oid.id = self.sip.get().id()
        return oid

    def call_method(self, method, **kwargs):
        cdef VariantMap parameters

        for name in kwargs:
            parameters[to_char_pointer(name)] = python_object_to_variant(
                kwargs[name])

        res = variant_to_python_object(
            self.sip.get().call_method(to_char_pointer(method), parameters))
        handle_errors("")
        return res

    def name(self):
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
    cdef Variant v
    cdef vector[Variant] vec
    cdef PObjectId oid

    if value is None:
        return Variant()

    # The order is important, the object character should
    # be preserved even if the PScriptInterface derived class
    # is iterable.
    if isinstance(value, PScriptInterface):
        # Map python object do id
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
        try:
            oid = get_value[ObjectId](value)
            ptr = get_instance(oid).lock()
            if ptr != shared_ptr[ScriptInterfaceBase]():
                so_name = to_str(ptr.get().name())
                if not so_name:
                    raise Exception(
                        "Script object without name returned from the core")
                # Fallback class, if nothing more specific is registered
                # for the script object name
                pclass = ScriptInterfaceHelper
                # Look up class
                if so_name in _python_class_by_so_name:
                    pclass = _python_class_by_so_name[so_name]
                pobj = pclass()
                poid = PObjectId()
                poid.id = ptr.get().id()
                pobj.set_sip_via_oid(poid)
                return pobj
            else:
                return None
        except:
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
    _so_name = None
    _so_bind_methods = ()
    _so_creation_policy = "GLOBAL"

    def __init__(self, **kwargs):
        super(ScriptInterfaceHelper, self).__init__(
            self._so_name, policy=self._so_creation_policy, **kwargs)
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

# Map from script object names to corresponding python classes
_python_class_by_so_name = {}


def script_interface_register(c):
    """Decorator used to register script interface classes
       This will store a name<->class relationship in a registry, so that parameters
       of type object can be instantiated as the correct python class
    """
    if not hasattr(c, "_so_name"):
        raise Exception(
            "Python classes representing a script object must define an _so_name attribute at class level")
    _python_class_by_so_name[c._so_name] = c
    return c


def init():
    initialize()
