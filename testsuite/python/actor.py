#
# Copyright (C) 2018-2019 The ESPResSo project
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


class TestActor(espressomd.actors.Actor):

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


class ActorTest(ut.TestCase):

    def test_ctor(self):
        a = TestActor(a=False, c=False)
        self.assertFalse(a.is_active())
        self.assertEqual(a.get_params(), a.default_params())
        self.assertEqual(a.system, None)

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
        actor = TestActor(a=False, c=False)
        self.actors.add(actor)
        with self.assertRaises(espressomd.highlander.ThereCanOnlyBeOne):
            self.actors.add(actor)
        # an actor can only be removed once
        self.actors.remove(actor)
        with self.assertRaises(Exception):
            self.actors.remove(actor)


if __name__ == "__main__":
    ut.main()
