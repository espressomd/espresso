include "myconfig.pxi"
from highlander import ThereCanOnlyBeOne


cdef class Actor:

    # Keys in active_list have to match the method name.
    active_list = dict(ElectrostaticInteraction=False,
                       MagnetostaticInteraction=False,
                       MagnetostaticExtension=False,
                       HydrodynamicInteraction=False,
                       ElectrostaticExtensions=False)

    # __getstate__ and __setstate__ define the pickle interaction
    def __getstate__(self):
        odict = self._params.copy()
        return odict

    def __setstate__(self,params):
        self._params=params
        self._set_params_in_es_core()

    def __init__(self, *args, **kwargs):
        self._isactive = False
        self._params = self.default_params()
        self.system = None

        # Check if all required keys are given
        for k in self.required_keys():
            if k not in kwargs:
                raise ValueError(
                    "At least the following keys have to be given as keyword arguments: " + self.required_keys().__str__() + " got " + kwargs.__str__())
            self._params[k] = kwargs[k]

        for k in kwargs:
            if k in self.valid_keys():
                self._params[k] = kwargs[k]
            else:
                raise KeyError("%s is not a vaild key" % k)
        self._set_params_in_es_core()

    def _activate(self):
        inter = self._get_interaction_type()
        if Actor.active_list[inter]:
            raise ThereCanOnlyBeOne(self.__class__.__bases__[0])
        Actor.active_list[inter] = True
        self.validate_params()
        self._activate_method()
        self._isactive = True

    def _deactivate(self):
        self._deactivate_method()
        self._isactive = False

    def is_valid(self):
        """Check, if the data stored in the instance still matches what is in Espresso"""
        # check, if the parameters saved in the class still match those
        # saved in Espresso
        temp_params = self._get_params_from_es_core()
        if self._params != temp_params:
            return False

        # If we're still here, the instance is valid
        return True

    def get_params(self):
        """Get interaction parameters"""
        # If this instance refers to an actual interaction defined in the es core, load
        # current parameters from there
        update = self._get_params_from_es_core()
        self._params.update(update)
        return self._params

    def set_params(self, **p):
        """Update parameters. Only given """
        # Check, if any key was passed, which is not known
        for k in p.keys():
            if k not in self.valid_keys():
                raise ValueError(
                    "Only the following keys are supported: " + self.valid_keys().__str__())

        # When an interaction is newly activated, all required keys must be
        # given
        if not self.is_active():
            for k in self.required_keys():
                if k not in p:
                    raise ValueError(
                        "At least the following keys have to be given as keyword arguments: " + self.required_keys().__str__())

        self._params.update(p)
        # vaidate updated parameters
        self.validate_params()
        # Put in values given by the user
        self._set_params_in_es_core()

    def _get_interaction_type(self):
        bases = self.class_lookup(self.__class__)
        for i in range(len(bases)):
            if bases[i].__name__ in Actor.active_list:
                return bases[i].__name__

    def class_lookup(self, cls):
        c = list(cls.__bases__)
        for base in c:
            c.extend(self.class_lookup(base))
        return c

    def is_active(self):
        print self.__class__.__name__, self._isactive
        return self._isactive

    def valid_keys(self):
        raise Exception(
            "Subclasses of %s must define the valid_keys() method." % self._get_interaction_type())

    def required_keys(self):
        raise Exception(
            "Subclasses of %s must define the required_keys() method." % self._get_interaction_type())

    def validate_params(self):
        raise Exception(
            "Subclasses of %s must define the validate_params() method." % self._get_interaction_type())

    def _get_params_from_es_core(self):
        raise Exception(
            "Subclasses of %s must define the _get_params_from_es_core() method." % self._get_interaction_type())

    def _set_params_in_es_core(self):
        raise Exception(
            "Subclasses of %s must define the _set_params_in_es_core() method." % self._get_interaction_type())

    def default_params(self):
        raise Exception(
            "Subclasses of %s must define the default_params() method." % self._get_interaction_type())

    def _activate_method(self):
        raise Exception(
            "Subclasses of %s must define the _activate_method() method." % self._get_interaction_type())

    def _deactivate_method(self):
        raise Exception(
            "Subclasses of %s must define the _deactivate_method() method." % self._get_interaction_type())


class Actors:

    active_actors = []

    def __init__(self, _system=None):
        self.system = _system

    def add(self, actor):
        if not actor in Actors.active_actors:
            actor.system = self.system
            Actors.active_actors.append(actor)
            actor._activate()
        else:
            raise ThereCanOnlyBeOne(actor)

    def __str__(self):
        print "Active Actors:"
        for actor in Actors.active_actors:
            print actor
        return ""

    def get(self):
        # for actor in Actors.activeActors:
        #     print actor.__class__.__name__
        return "%s" % Actors.active_actors

    def deactivate(self, actor):
        actor._deactivate()
