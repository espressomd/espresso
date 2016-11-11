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

    cdef Variant make_variant(self, ParameterType type, value):
        if < int > type == <int > BOOL:
            return make_variant[bool]( < bool > value)
        if < int > type == <int > INT:
            return make_variant[int]( < int > value)
        if < int > type == <int > DOUBLE:
            return make_variant[double]( < double > value)
        if < int > type == <int > STRING:
            return make_variant[string](to_char_pointer(value))
        if < int > type == <int > INT_VECTOR:
            return make_variant[vector[int] ]( < vector[int] > value)
        if < int > type == <int > DOUBLE_VECTOR:
            return make_variant[vector[double] ]( < vector[double] > value)
        if < int > type == <int > VECTOR3D:
            return make_variant[Vector3d](Vector3d( < vector[double] > value))
        if < int > type == <int > VECTOR2D:
            return make_variant[Vector2d](Vector2d( < vector[double] > value))

        raise Exception("Unkown type")

    def id(self):
        oid = PObjectId()
        oid.id = self.sip.get().id()
        return oid

    def call_method(self, method, **kwargs):
        cdef map[string, Variant] parameters
        cdef PObjectId oid

        for name in kwargs:
            if isinstance(kwargs[name], PScriptInterface):
                # Map python object do id
                oid = kwargs[name].id()
                parameters[to_char_pointer(name)] = oid.id
            else:
                parameters[to_char_pointer(name)] = kwargs[name]

        return self.variant_to_python_object(self.sip.get().call_method(to_char_pointer(method), parameters))

    def set_params(self, **kwargs):
        cdef ParameterType type
        cdef map[string, Variant] parameters
        cdef PObjectId oid
        cdef string name

        for pname in kwargs:
            name = to_char_pointer(pname)

            try:
                type = self.parameters.at(name).type()
            except:
                raise ValueError("Unknown parameter %s" % name)

            # Check number of elements if applicable
            if < int > type in [ < int > INT_VECTOR, < int > DOUBLE_VECTOR, < int > VECTOR2D, < int > VECTOR3D]:
                n_elements = self.parameters[name].n_elements()
                if n_elements!=0 and not (len(kwargs[pname]) == n_elements):
                    raise ValueError(
                        "Value of %s expected to be %i elements" % (name, n_elements))

            # Objects have to be translated to ids
            if < int > type is < int > OBJECT:
                oid = kwargs[pname].id()
                parameters[name] = oid.id
            else:
                parameters[name] = self.make_variant(type, kwargs[pname])

        self.sip.get().set_parameters(parameters)

    cdef variant_to_python_object(self, Variant value):
        cdef ObjectId oid
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
            try:
                oid = get[ObjectId](value)
                pobj = PScriptInterface()
                pobj.set_sip(get_instance(oid).lock())
                return pobj
            except:
                return None

        raise Exception("Unkown type")

    def get_parameter(self, name):
        cdef ParameterType type = self.parameters[name].type()
        cdef Variant value = self.sip.get().get_parameter(name)
        return self.variant_to_python_object(value)

    def get_params(self):
        odict = {}
        for pair in self.parameters:
            try:
                odict[pair.first] = self.get_parameter(pair.first)
            except:
                pass
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
