cdef class PScriptInterface:
    def __init__(self, name=None):
        if name:
            self.sip = make_shared(name)
            self.parameters = self.sip.get().all_parameters()
        else:
            self.sip = shared_ptr[ScriptInterfaceBase]()
            self.parameters = map[string, Parameter]()

    def __richcmp__(a, b, op):
        if op == 2:
            return a.id() == b.id()
        else:
            raise NotImplementedError

    def id(self):
        return self.sip.get().id()

    def _ref_count(self):
        return self.sip.use_count()

    cdef set_sip(self, shared_ptr[ScriptInterfaceBase] sip):
        self.sip = sip
        self.parameters = self.sip.get().all_parameters()

    cdef Variant make_variant(self, ParameterType type, value):
        if < int > type == <int > BOOL:
            return make_variant[bool]( < bool > value)
        if < int > type == <int > INT:
            return make_variant[int]( < int > value)
        if < int > type == <int > DOUBLE:
            return make_variant[double]( < double > value)
        if < int > type == <int > STRING:
            return make_variant[string]( < string > value)
        if < int > type == <int > INT_VECTOR:
            return make_variant[vector[int] ]( < vector[int] > value)
        if < int > type == <int > DOUBLE_VECTOR:
            return make_variant[vector[double] ]( < vector[double] > value)
        if < int > type == <int > VECTOR3D:
            return make_variant[Vector3d](Vector3d( < vector[double] > value))
        if < int > type == <int > VECTOR2D:
            return make_variant[Vector2d](Vector2d( < vector[double] > value))

        raise Exception("Unkown type")

    def call_method(self, method, **kwargs):
        cdef map[string, Variant] parameters

        for name in kwargs:
            if isinstance(kwargs[name], PScriptInterface):
                # Map python object do id
                parameters[name] = OId(kwargs[name].id())
            else:
                parameters[name] = kwargs[name]

        return self.variant_to_python_object(self.sip.get().call_method(method, parameters))


    def set_parameters(self, **kwargs):
        cdef ParameterType type
        cdef map[string, Variant] parameters

        for name in kwargs:
            try:
                self.parameters.at(name)
            except:
                raise ValueError("Unknown parameter %s" % name)

            type = self.parameters[name].type()

            # Check number of elements if applicable
            if < int > type in [ < int > INT_VECTOR, < int > DOUBLE_VECTOR, < int > VECTOR2D, < int > VECTOR3D]:
                n_elements = self.parameters[name].n_elements()
                if not (len(kwargs[name]) == n_elements):
                    raise ValueError(
                        "Value of %s expected to be %i elements" % (name, n_elements))

            # Objects have to be translated to ids
            if < int > type is < int > OBJECT:
                parameters[name] = OId(kwargs[name].id())
            else:
                parameters[name] = self.make_variant(type, kwargs[name])

        self.sip.get().set_parameters(parameters)

    cdef variant_to_python_object(self, Variant value):
        cdef int type = value.which()
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
        if < int > type == <int > VECTOR3D:
            return get[Vector3d](value).as_vector()
        if < int > type == <int > VECTOR2D:
            return get[Vector2d](value).as_vector()
        if < int > type == <int > OBJECT:
            # Get the id and build a curresponding object
            val = get[OId](value).id
            if val >= 0:
                pobj = PScriptInterface()
                pobj.set_sip(get_instance(val).lock())
                return pobj
            else:
                return None

        raise Exception("Unkown type")

    def get_parameter(self, name):
        cdef ParameterType type = self.parameters[name].type()
        cdef Variant value = self.sip.get().get_parameter(name)
        return self.variant_to_python_object(value)

    def get_parameters(self):
        odict = {}
        for pair in self.parameters:
            odict[pair.first] = self.get_parameter(pair.first)
        return odict

initialize()
