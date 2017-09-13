from espressomd.utils import to_char_pointer,to_str

cdef class PObjectId(object):
    cdef ObjectId id

    def __richcmp__(PObjectId a, PObjectId b, op):
        if op == 2:
            return a.id == b.id
        else:
            raise NotImplementedError

cdef class PScriptInterface(object):
    def __init__(self, name=None, policy="GLOBAL"):
        if name:
            if(policy=="GLOBAL"):
                self.sip = make_shared(to_char_pointer(name), GLOBAL)
            else:
                self.sip = make_shared(to_char_pointer(name), LOCAL)
            self.parameters = self.sip.get().valid_parameters()
        else:
            raise Exception("the name parameter has to be set.")

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

    def set_sip_via_oid(self,PObjectId id):
        """Set the shared_ptr to the script object in the core via the object id"""
        oid=id.id
        try:
            ptr = get_instance(oid).lock()
            self.set_sip(ptr)
        except:
            raise Exception("Could not get sip for given_id")
    
    def _valid_parameters(self):
        parameters = []

        for p in self.sip.get().valid_parameters():
            parameters.append(to_str(p.first))

        return parameters

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
        return to_str(self.sip.get().name())

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

    cdef variant_to_python_object(self, Variant value) except +:
        cdef ObjectId oid
        cdef vector[Variant] vec
        cdef int type = value.which()
        cdef shared_ptr[ScriptInterfaceBase] ptr

        if < int > type == <int > NONE:
            return None
        if < int > type == <int > BOOL:
            return get[bool](value)
        if < int > type == <int > INT:
            return get[int](value)
        if < int > type == <int > DOUBLE:
            return get[double](value)
        if < int > type == <int > STRING:
            return to_str(get[string](value))
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
                    so_name=to_str(ptr.get().name())
                    # Fallback class, if nothing omre specific is registered for the script object name
                    pclass=ScriptInterfaceHelper
                    # Look up class
                    if so_name in _python_class_by_so_name:
                        pclass =_python_class_by_so_name[so_name]
                    pobj = pclass()
                    poid=PObjectId()
                    poid.id=ptr.get().id()
                    pobj.set_sip_via_oid(poid)
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
        cdef Variant value = self.sip.get().get_parameter(to_char_pointer(name))
        return self.variant_to_python_object(value)

    def get_params(self):
        cdef map[string, Variant] params = self.sip.get().get_parameters()
        odict = {}
        for pair in params:
            odict[to_str(pair.first)] = self.variant_to_python_object(pair.second)

        return odict

class ScriptInterfaceHelper(PScriptInterface):
    _so_name = None
    _so_bind_methods =()
    _so_creation_policy = "GLOBAL"

    def __init__(self, **kwargs):
        super(ScriptInterfaceHelper,self).__init__(self._so_name, self._so_creation_policy)
        self.set_params(**kwargs)
        self.define_bound_methods()

    def __dir__(self):
        return self.__dict__.keys() + self._valid_parameters()

    def __getattr__(self, attr):
        if attr in self._valid_parameters():
            return self.get_parameter(to_char_pointer(attr))
        else:
            try:
                return self.__dict__[attr]
            except KeyError:
                raise AttributeError

    def __setattr__(self, attr, value):
        if attr in self._valid_parameters():
            self.set_params(**{attr:value})
        else:
            self.__dict__[attr] = value

    def generate_caller(self,method_name):
        def template_method(**kwargs):
            res=self.call_method(method_name,**kwargs)
            return res

        return template_method

    def define_bound_methods(self):
        for method_name in self._so_bind_methods:
            setattr(self,method_name,self.generate_caller(method_name))





# Map from script object names to corresponding python classes
_python_class_by_so_name ={}


def script_interface_register(c):
    """Decorator used to register script interface classes
       This will store a name<->class relationship in a registry, so that parameters
       of type object can be instanciated as the correct python class
    """
    if not hasattr(c,"_so_name"):
        raise Exception("Python classes representing a script object must define an _so_name attribute at class level")
    _python_class_by_so_name[c._so_name]=c
    return c
    

def init():
    initialize()
