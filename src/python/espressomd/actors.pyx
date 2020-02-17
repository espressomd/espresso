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
include "myconfig.pxi"
from .highlander import ThereCanOnlyBeOne
from .utils cimport handle_errors

cdef class Actor:

    """
    Abstract base class for interactions affecting particles in the system,
    such as electrostatics, magnetostatics or LB fluids. Derived classes must
    implement the interface to the relevant core objects and global variables.
    """

    # Keys in active_list have to match the method name.
    active_list = dict(ElectrostaticInteraction=False,
                       MagnetostaticInteraction=False,
                       MagnetostaticExtension=False,
                       HydrodynamicInteraction=False,
                       Scafacos=False)

    # __getstate__ and __setstate__ define the pickle interaction
    def __getstate__(self):
        odict = self._params.copy()
        return odict

    def __setstate__(self, params):
        self._params = params
        self._set_params_in_es_core()

    def __init__(self, *args, **kwargs):
        self._isactive = False
        self._params = self.default_params()

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
                raise KeyError("%s is not a valid key" % k)

    def _activate(self):
        inter = self._get_interaction_type()
        if inter in Actor.active_list:
            if Actor.active_list[inter]:
                raise ThereCanOnlyBeOne(self.__class__.__bases__[0])
            Actor.active_list[inter] = True

        self.validate_params()
        self._activate_method()
        handle_errors("Activation of an actor")
        self._isactive = True

    def _deactivate(self):
        self._deactivate_method()
        handle_errors("Deactivation of an actor")
        self._isactive = False
        inter = self._get_interaction_type()
        if inter in Actor.active_list:
            if not Actor.active_list[inter]:
                raise Exception(
                    "Class not registered in Actor.active_list " + self.__class__.__bases__[0])
            Actor.active_list[inter] = False

    def is_valid(self):
        """
        Check if the data stored in this instance still matches the
        corresponding data in the core.
        """
        temp_params = self._get_params_from_es_core()
        if self._params != temp_params:
            return False

        # If we're still here, the instance is valid
        return True

    def get_params(self):
        """Get interaction parameters"""
        # If this instance refers to an actual interaction defined in the es
        # core, load current parameters from there
        if self.is_active():
            update = self._get_params_from_es_core()
            self._params.update(update)
        return self._params

    def set_params(self, **p):
        """Update the given parameters."""
        # Check if keys are valid
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
        # validate updated parameters
        self.validate_params()
        # Put in values given by the user
        if self.is_active():
            self._set_params_in_es_core()

    def __str__(self):
        return self.__class__.__name__ + "(" + str(self.get_params()) + ")"

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
        return self._isactive

    def valid_keys(self):
        """Virtual method."""
        raise Exception(
            "Subclasses of %s must define the valid_keys() method." % self._get_interaction_type())

    def required_keys(self):
        """Virtual method."""
        raise Exception(
            "Subclasses of %s must define the required_keys() method." % self._get_interaction_type())

    def validate_params(self):
        """Virtual method."""
        raise Exception(
            "Subclasses of %s must define the validate_params() method." % self._get_interaction_type())

    def _get_params_from_es_core(self):
        """Virtual method."""
        raise Exception(
            "Subclasses of %s must define the _get_params_from_es_core() method." % self._get_interaction_type())

    def _set_params_in_es_core(self):
        """Virtual method."""
        raise Exception(
            "Subclasses of %s must define the _set_params_in_es_core() method." % self._get_interaction_type())

    def default_params(self):
        """Virtual method."""
        raise Exception(
            "Subclasses of %s must define the default_params() method." % self._get_interaction_type())

    def _activate_method(self):
        """Virtual method."""
        raise Exception(
            "Subclasses of %s must define the _activate_method() method." % self._get_interaction_type())

    def _deactivate_method(self):
        """Virtual method."""
        raise Exception(
            "Subclasses of %s must define the _deactivate_method() method." % self._get_interaction_type())


class Actors:

    """
    Container for :class:`Actor` objects.
    """

    active_actors = []

    def __getstate__(self):
        return self.active_actors

    def __setstate__(self, active_actors):
        self.active_actors[:] = []
        for a in active_actors:
            self.active_actors.append(a)
            a._activate()

    def add(self, actor):
        """
        Parameters
        ----------
        actor : :class:`Actor`
            Actor to add to this container.

        """
        if actor not in Actors.active_actors:
            self.active_actors.append(actor)
            actor._activate()
        else:
            raise ThereCanOnlyBeOne(actor)

    def remove(self, actor):
        """
        Parameters
        ----------
        actor : :class:`Actor`
            Actor to remove from this container.

        """
        if actor not in self.active_actors:
            raise Exception("Actor is not active")
        actor._deactivate()
        self.active_actors.remove(actor)

    def clear(self):
        """Remove all actors."""
        for a in self.active_actors:
            self.remove(a)

    def __str__(self):
        return str(self.active_actors)

    def __getitem__(self, key):
        return self.active_actors[key]

    def __len__(self):
        return len(self.active_actors)

    def __iter__(self):
        for a in self.active_actors:
            yield a

    def __delitem__(self, idx):
        actor = self[idx]
        self.remove(actor)
