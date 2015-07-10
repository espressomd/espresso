include "myconfig.pxi"
from highlander import ThereCanOnlyBeOne


cdef class Actor:
    activeList = dict(ElectrostaticInteraction=False,
                      MagnetostaticInteraction=False,
                      MagnetostaticExtension=False,
                      HydrodynamicInteraction=False,
                      ElectrostaticExtensions=False)

    def __cinit__(self, *args, **kwargs):
        self._isactive = False
        self._params = self.defaultParams()
        self.system = None

        # Check if all required keys are given
        for k in self.requiredKeys():
            if k not in kwargs:
                raise ValueError(
                    "At least the following keys have to be given as keyword arguments: " + self.requiredKeys().__str__() + " got " + kwargs.__str__())
            self._params[k] = kwargs[k]

        for k in kwargs:
            if k in self.validKeys():
                self._params[k] = kwargs[k]
            else:
                raise KeyError("%s is not a vaild key" % k)

    def _activate(self):
        inter = self._getInteractionType()
        if Actor.activeList[inter]:
            raise ThereCanOnlyBeOne(self.__class__.__bases__[0])
        Actor.activeList[inter] = True
        self.validateParams()
        self._activateMethod()
        self._isactive = True

    def _deactivate(self):
        self._deactivateMethod()
        self._isactive = False

    def isValid(self):
        """Check, if the data stored in the instance still matches what is in Espresso"""
        # check, if the bond parameters saved in the class still match those
        # saved in Espresso
        tempParams = self._getParamsFromEsCore()
        if self._params != tempParams:
            return False

        # If we're still here, the instance is valid
        return True

    def getParams(self):
        """Get interaction parameters"""
        # If this instance refers to an actual interaction defined in the es core, load
        # current parameters from there
        update = self._getParamsFromEsCore()
        self._params.update(update)
        return self._params

    def setParams(self, **p):
        """Update parameters. Only given """
        # Check, if any key was passed, which is not known
        for k in p.keys():
            if k not in self.validKeys():
                raise ValueError(
                    "Only the following keys are supported: " + self.validKeys().__str__())

        # When an interaction is newly activated, all required keys must be
        # given
        if not self.isActive():
            for k in self.requiredKeys():
                if k not in p:
                    raise ValueError(
                        "At least the following keys have to be given as keyword arguments: " + self.requiredKeys().__str__())

        self._params.update(p)
        # vaidate updated parameters
        self.validateParams()
        # Put in values given by the user
        self._setParamsInEsCore()

    def _getInteractionType(self):
        return self.__class__.__bases__[0].__name__

    def isActive(self):
        print self.__class__.__name__, self._isactive
        return self._isactive

    def validKeys(self):
        raise Exception(
            "Subclasses of %s must define the validKeys() method." % self._getInteractionType())

    def requiredKeys(self):
        raise Exception(
            "Subclasses of %s must define the requiredKeys() method." % self._getInteractionType())

    def validateParams(self):
        raise Exception(
            "Subclasses of %s must define the validateParams() method." % self._getInteractionType())

    def _getParamsFromEsCore(self):
        raise Exception(
            "Subclasses of %s must define the _getParamsFromEsCore() method." % self._getInteractionType())

    def _setParamsInEsCore(self):
        raise Exception(
            "Subclasses of %s must define the _setParamsInEsCore() method." % self._getInteractionType())

    def defaultParams(self):
        raise Exception(
            "Subclasses of %s must define the defaultParams() method." % self._getInteractionType())

    def _activateMethod(self):
        raise Exception(
            "Subclasses of %s must define the _activateMethod() method." % self._getInteractionType())

    def _deactivateMethod(self):
        raise Exception(
            "Subclasses of %s must define the _deactivateMethod() method." % self._getInteractionType())


class Actors:
    activeActors = []

    def __init__(self, _system=None):
        self.system = _system

    def add(self, actor):
        if not actor in Actors.activeActors:
            actor.system = self.system
            Actors.activeActors.append(actor)
            actor._activate()
        else:
            raise ThereCanOnlyBeOne(actor)

    def __str__(self):
        print "Active Actors:"
        for actor in Actors.activeActors:
            print actor
        return ""

    def get(self):
        # for actor in Actors.activeActors:
        #     print actor.__class__.__name__
        return "%s" % Actors.activeActors

    def deactivate(self, actor):
        actor._deactivate()
