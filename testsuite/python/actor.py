#
# Copyright (C) 2018-2022 The ESPResSo project
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
#

"""
Testmodule for the actor base class.

"""

import unittest as ut
import espressomd.actors
import espressomd.highlander
import espressomd.utils as utils


class BaseActor:

    """
    Abstract base class for interactions affecting particles in the system,
    such as LB fluids. Derived classes must implement the interface to the
    relevant core objects and global variables.
    """

    # Keys in active_list have to match the method name.
    active_list = dict(HydrodynamicInteraction=False)

    def __init__(self, **kwargs):
        self._isactive = False
        utils.check_valid_keys(self.valid_keys(), kwargs.keys())
        utils.check_required_keys(self.required_keys(), kwargs.keys())
        self._params = self.default_params()
        self._params.update(kwargs)

    def _activate(self):
        inter = self._get_interaction_type()
        if inter in BaseActor.active_list:
            if BaseActor.active_list[inter]:
                raise espressomd.highlander.ThereCanOnlyBeOne(
                    self.__class__.__bases__[0])
            BaseActor.active_list[inter] = True

        self.validate_params()
        self._activate_method()
        utils.handle_errors("Activation of an actor")
        self._isactive = True

    def _deactivate(self):
        self._deactivate_method()
        utils.handle_errors("Deactivation of an actor")
        self._isactive = False
        inter = self._get_interaction_type()
        if inter in BaseActor.active_list:
            if not BaseActor.active_list[inter]:
                raise Exception(
                    f"Class not registered in Actor.active_list: {self.__class__.__bases__[0].__name__}")
            BaseActor.active_list[inter] = False

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
        utils.check_valid_keys(self.valid_keys(), p.keys())

        # When an interaction is newly activated, all required keys must be
        # given
        if not self.is_active():
            utils.check_required_keys(self.required_keys(), p.keys())

        self._params.update(p)
        # validate updated parameters
        self.validate_params()
        # Put in values given by the user
        if self.is_active():
            self._set_params_in_es_core()

    def is_active(self):
        return self._isactive


class TestActor(BaseActor):

    def __init__(self, *args, **kwargs):
        self._core_args = None
        self._activated = False
        self._deactivated = False
        self._validated = False

        super().__init__(*args, **kwargs)

    def _get_params_from_es_core(self):
        return self._core_args

    def _set_params_in_es_core(self):
        self._core_args = self._params

    def valid_keys(self):
        return {"a", "b", "c"}

    def required_keys(self):
        return {"a", "c"}

    def default_params(self):
        return {"a": False, "b": False, "c": False}

    def _activate_method(self):
        self._activated = True

    def _deactivate_method(self):
        self._deactivated = True

    def validate_params(self):
        self._validated = True

    def _get_interaction_type(self):
        return None


class TestHydrodynamicActor(TestActor):

    def _get_interaction_type(self):
        return "HydrodynamicInteraction"


class ActorTest(ut.TestCase):

    def test_ctor(self):
        a = TestActor(a=False, c=False)
        self.assertFalse(a.is_active())
        self.assertEqual(a.get_params(), a.default_params())

    def test_params_non_active(self):
        a = TestActor(a=True, c=True)

        a.set_params(a=False, b=True, c=False)
        params = a.get_params()
        self.assertEqual(params["a"], False)
        self.assertEqual(params["b"], True)
        self.assertEqual(params["c"], False)
        self.assertEqual(a._core_args, None)

    def test_params_active(self):
        a = TestActor(a=True, c=True)
        a._activate()

        a.set_params(a=False, b=True, c=False)
        params = a.get_params()
        self.assertEqual(params["a"], False)
        self.assertEqual(params["b"], True)
        self.assertEqual(params["c"], False)
        self.assertEqual(a._core_args, params)

    def test_activation(self):
        a = TestActor(a=True, c=True)
        a._activate()
        self.assertTrue(a.is_active())

    def test_deactivation(self):
        a = TestActor(a=True, c=True)
        a._activate()
        self.assertTrue(a.is_active())
        a._deactivate()
        self.assertFalse(a.is_active())

        params = a.get_params()
        self.assertEqual(params["a"], True)
        self.assertEqual(params["b"], False)
        self.assertEqual(params["c"], True)

    def test_exception(self):
        error_msg_valid = (r"Only the following keys can be given as keyword arguments: "
                           r"\['a', 'b', 'c'\], got \['a', 'c', 'd'\] \(unknown \['d'\]\)")
        error_msg_required = (r"The following keys have to be given as keyword arguments: "
                              r"\['a', 'c'\], got \['a'\] \(missing \['c'\]\)")
        with self.assertRaisesRegex(ValueError, error_msg_valid):
            TestActor(a=True, c=True, d=True)
        with self.assertRaisesRegex(ValueError, error_msg_required):
            TestActor(a=True)
        valid_actor = TestActor(a=True, c=True)
        with self.assertRaisesRegex(ValueError, error_msg_valid):
            valid_actor.set_params(a=True, c=True, d=True)
        with self.assertRaisesRegex(ValueError, error_msg_required):
            valid_actor.set_params(a=True)


class ActorsTest(ut.TestCase):

    actors = espressomd.actors.Actors()

    def tearDown(self):
        self.actors.clear()

    def test_clear(self):
        # clearing the list of actors removes all of them
        for actors_size in range(10):
            for _ in range(actors_size):
                actor = TestActor(a=False, c=False)
                self.actors.add(actor)
            self.assertEqual(len(self.actors), actors_size)
            self.actors.clear()
            self.assertEqual(len(self.actors), 0)

    def test_deactivation(self):
        actor = TestActor(a=False, c=False)
        self.assertFalse(actor.is_active())
        # adding an actor activates it
        self.actors.add(actor)
        self.assertTrue(actor.is_active())
        # removing an actor deactivates it
        self.actors.clear()
        self.assertFalse(actor.is_active())
        # re-adding an actor re-activates it
        self.actors.add(actor)
        self.assertTrue(actor.is_active())
        # removing an actor deactivates it
        del self.actors[0]
        self.assertFalse(actor.is_active())

    def test_unique(self):
        # an actor can only be added once
        actor = TestHydrodynamicActor(a=False, c=False)
        self.actors.add(actor)
        with self.assertRaises(espressomd.highlander.ThereCanOnlyBeOne):
            self.actors.add(actor)
        with self.assertRaises(espressomd.highlander.ThereCanOnlyBeOne):
            actor._activate()
        # an actor can only be removed once
        self.actors.remove(actor)
        with self.assertRaisesRegex(Exception, "Actor is not active"):
            self.actors.remove(actor)
        with self.assertRaisesRegex(Exception, "Class not registered.*: TestActor"):
            actor._deactivate()


if __name__ == "__main__":
    ut.main()
