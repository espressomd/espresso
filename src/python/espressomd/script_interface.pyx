from espressomd.utils import to_char_pointer

cdef class PScriptInterface:
    def __init__(self, name=None):
        if name:
            self.sip = make_shared(to_char_pointer(name))
            self.parameters = self.sip.get().valid_parameters()
        else:
            self.sip = shared_ptr[ScriptInterfaceBase]()
            self.parameters = map[string, Parameter]()

    def __richcmp__(a, b, op):
        if op == 2:
            return a.id() == b.id()
        else:
            raise NotImplementedError

    def _ref_count(self):
        return self.sip.use_count()

    cdef set_sip(self, shared_ptr[ScriptInterfaceBase] sip):
        self.sip = sip
        self.parameters = self.sip.get().valid_parameters()

    def id(self):
        oid = PObjectId()
        oid.id = self.sip.get().id()
        return oid

    def call_method(self, method, **kwargs):
        cdef map[string, Variant] parameters

        for name in kwargs:
            parameters[to_char_pointer(name)] = self.python_object_to_variant(kwargs[name])

        return self.variant_to_python_object(self.sip.get().call_method(to_char_pointer(method), parameters))

    def name(self):
        return self.sip.get().name()

    def set_params(self, **kwargs):
        cdef ParameterType type
        cdef map[string, Variant] parameters
        cdef Variant v

        for pname in kwargs:
            name = to_char_pointer(pname)

            try:
                type = self.parameters.at(name).type()
            except:
                raise ValueError("Unknown parameter %s" % name)

            # Check number of elements if applicable
            if < int > type in [ < int > INT_VECTOR, < int > DOUBLE_VECTOR]:
                n_elements = self.parameters[name].n_elements()
                if n_elements!=0 and not (len(kwargs[pname]) == n_elements):
                    raise ValueError(
                        "Value of %s expected to be %i elements" % (name, n_elements))

            # We accept floats for ints (but not the other way round)
            if <int> type == <int> DOUBLE and isinstance(kwargs[pname], int):
                kwargs[pname] = float(kwargs[pname])
            # We already know that the argument is an iterable of the correct length
            elif <int> type == <int> DOUBLE_VECTOR:
                for i in range(len(kwargs[pname])):
                    if isinstance(kwargs[pname][i], int):
                        kwargs[pname][i] = float(kwargs[pname][i])

            v = self.python_object_to_variant(kwargs[pname])

            if v.which() == <int> type:
                parameters[name] = v
            else:
                raise ValueError("Wrong type for parameter '%s': Expected %s, but got %s" % (pname, get_type_label(type), get_type_label(v)))

        self.sip.get().set_parameters(parameters)

    cdef Variant python_object_to_variant(self, value):
        cdef Variant v
        cdef vector[Variant] vec
        cdef PObjectId oid

        # The order is important, the object character should
        # be preserved even if the PScriptInterface derived class
        # is iterable.
        if isinstance(value, PScriptInterface):
            # Map python object do id
            oid = value.id()
            return make_variant[ObjectId](oid.id)
        elif hasattr(value, '__iter__') and not(type(value) == str):
            for e in value:
                vec.push_back(self.python_object_to_variant(e))
            v = make_variant[vector[Variant]](vec)
            transform_vectors(v)
            return v
        elif type(value) == str:
            return make_variant[string](to_char_pointer(value))
        elif type(value) == int:
            return make_variant[int](value)
        elif type(value) == float:
            return make_variant[float](value)
        elif type(value) == type(True):
            return make_variant[bool](value)
        else:
            raise TypeError("Unkown type for conversion to Variant")

    cdef variant_to_python_object(self, Variant value):
        cdef ObjectId oid
        cdef vector[Variant] vec
        cdef int type = value.which()
        cdef shared_ptr[ScriptInterfaceBase] ptr

        if < int > type == <int > BOOL:
            return get[bool](value)
        if < int > type == <int > INT:
            return get[int](value)
        if < int > type == <int > DOUBLE:
            return get[double](value)
        if < int > type == <int > STRING:
            return get[string](value)
        if < int > type == <int > INT_VECTOR:
            return get[vector[int]](value)
        if < int > type == <int > DOUBLE_VECTOR:
            return get[vector[double]](value)
        if < int > type == <int > OBJECTID:
            # Get the id and build a curresponding object
            try:
                oid = get[ObjectId](value)
                ptr = get_instance(oid).lock()
                if ptr != shared_ptr[ScriptInterfaceBase]():
                    pobj = PScriptInterface()
                    pobj.set_sip(ptr)
                    return pobj
                else:
                    return None
            except:
                return None
        if < int > type == < int > VECTOR:
            vec = get[vector[Variant]](value)
            res = []

            for i in vec:
                res.append(self.variant_to_python_object(i))

            return res

        raise Exception("Unkown type")

    def get_parameter(self, name):
        cdef Variant value = self.sip.get().get_parameter(name)
        return self.variant_to_python_object(value)

    def get_params(self):
        cdef map[string, Variant] params = self.sip.get().get_parameters()
        odict = {}
        for pair in params:
            odict[pair.first] = self.variant_to_python_object(pair.second)

        return odict

class ScriptInterfaceHelper(PScriptInterface):
    _so_name = None
    _so_bind_methods =()

    def __init__(self,**kwargs):
        super(ScriptInterfaceHelper,self).__init__(self._so_name)
        self.set_params(**kwargs)
        self.define_bound_methods()

    def generate_caller(self,method_name):
        def template_method(**kwargs):
            res=self.call_method(method_name,**kwargs)
            return res

        return template_method

    def define_bound_methods(self):
        for method_name in self._so_bind_methods:
            setattr(self,method_name,self.generate_caller(method_name))




initialize()
