#
# Copyright (C) 2013-2022 The ESPResSo project
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
import unittest as ut
import unittest_decorators as utx
import tests_common

import espressomd
import espressomd.interactions


class Test(ut.TestCase):
    system = espressomd.System(box_l=[30.0, 30.0, 30.0])

    def tearDown(self):
        if espressomd.has_features(["LENNARD_JONES"]):
            self.system.non_bonded_inter[0, 0].lennard_jones.set_params(
                epsilon=0., sigma=0., cutoff=0., shift=0.)

    def intersMatch(self, inType, outInter, inParams, outParams, msg_long):
        """Check, if the interaction type set and gotten back as well as the
        bond parameters set and gotten back match. Only check keys present in
        ``inParams``.
        """
        self.assertIsInstance(outInter, inType)
        tests_common.assert_params_match(self, inParams, outParams, msg_long)

    def parameterKeys(self, interObject):
        """
        Check :meth:`~espressomd.interactions.NonBondedInteraction.valid_keys`
        and :meth:`~espressomd.interactions.NonBondedInteraction.required_keys`
        return sets, and that
        :meth:`~espressomd.interactions.NonBondedInteraction.default_params`
        returns a dictionary with the correct keys.

        Parameters
        ----------
        interObject: instance of a class derived from :class:`espressomd.interactions.NonBondedInteraction`
            Object of the interaction to test, e.g.
            :class:`~espressomd.interactions.LennardJonesInteraction`
        """
        classname = interObject.__class__.__name__
        valid_keys = interObject.valid_keys()
        required_keys = interObject.required_keys()
        default_keys = set(interObject.default_params().keys())
        self.assertIsInstance(valid_keys, set,
                              "{}.valid_keys() must return a set".format(
                                  classname))
        self.assertIsInstance(required_keys, set,
                              "{}.required_keys() must return a set".format(
                                  classname))
        self.assertTrue(default_keys.issubset(valid_keys),
                        "{}.default_params() has unknown parameters: {}".format(
            classname, default_keys.difference(valid_keys)))
        self.assertTrue(default_keys.isdisjoint(required_keys),
                        "{}.default_params() has extra parameters: {}".format(
            classname, default_keys.intersection(required_keys)))
        self.assertSetEqual(default_keys, valid_keys - required_keys,
                            "{}.default_params() should have keys: {}, got: {}".format(
                                classname, valid_keys - required_keys, default_keys))

    def generateTestForNon_bonded_interaction(
            _partType1, _partType2, _interClass, _params, _interName):
        """Generates test cases for checking interaction parameters set and
        gotten back from the espresso core actually match those in the Python
        classes. Only keys which are present in ``_params`` are checked.

        Parameters
        ----------
        _partType1, _partType2: :obj:`int`
            Particle type ids to check on
        _interClass: class derived from :class:`espressomd.interactions.NonBondedInteraction`
            Class of the interaction to test, e.g.
            :class:`~espressomd.interactions.LennardJonesInteraction`
        _params: :obj:`dict`
            Interaction parameters, e.g. ``{"k": 1., "r_0": 0}``
        _interName: :obj:`str`
            Name of the interaction property to set (e.g. ``"lennard_jones"``)
        """
        partType1 = _partType1
        partType2 = _partType2
        interClass = _interClass
        params = _params
        interName = _interName

        def func(self):
            # This code is run at the execution of the generated function.
            # It will use the state of the variables in the outer function,
            # which was there, when the outer function was called

            # Set parameters
            getattr(self.system.non_bonded_inter[partType1, partType2],
                    interName).set_params(**params)

            # Read them out again
            outInter = getattr(
                self.system.non_bonded_inter[partType1, partType2], interName)
            outParams = outInter.get_params()

            self.intersMatch(
                interClass, outInter, params, outParams,
                "{}: value set and value gotten back differ for particle types {} and {}: {} vs. {}"
                .format(interClass(**params).type_name(), partType1, partType2,
                        params, outParams))
            self.parameterKeys(outInter)

        return func

    if espressomd.has_features(["LENNARD_JONES"]):
        test_lj1 = generateTestForNon_bonded_interaction(
            0, 0, espressomd.interactions.LennardJonesInteraction,
            {"epsilon": 1., "sigma": 2., "cutoff": 3.,
             "shift": 4., "offset": 5., "min": 7.},
            "lennard_jones")
        test_lj2 = generateTestForNon_bonded_interaction(
            0, 0, espressomd.interactions.LennardJonesInteraction,
            {"epsilon": 1.3, "sigma": 2.2, "cutoff": 3.4,
             "shift": 4.1, "offset": 5.1, "min": 7.1},
            "lennard_jones")
        test_lj3 = generateTestForNon_bonded_interaction(
            0, 0, espressomd.interactions.LennardJonesInteraction,
            {"epsilon": 1.3, "sigma": 2.2, "cutoff": 3.4,
             "shift": 4.1, "offset": 5.1, "min": 7.1},
            "lennard_jones")

    if espressomd.has_features(["LENNARD_JONES_GENERIC"]):
        test_ljgen1 = generateTestForNon_bonded_interaction(
            0, 0, espressomd.interactions.GenericLennardJonesInteraction,
            {"epsilon": 1., "sigma": 2., "cutoff": 3., "shift": 4.,
             "offset": 5., "e1": 7, "e2": 8, "b1": 9., "b2": 10.},
            "generic_lennard_jones")
        test_ljgen2 = generateTestForNon_bonded_interaction(
            0, 0, espressomd.interactions.GenericLennardJonesInteraction,
            {"epsilon": 1.1, "sigma": 2.1, "cutoff": 3.1, "shift": 4.1,
             "offset": 5.1, "e1": 71, "e2": 81, "b1": 9.1, "b2": 10.1},
            "generic_lennard_jones")
        test_ljgen3 = generateTestForNon_bonded_interaction(
            0, 0, espressomd.interactions.GenericLennardJonesInteraction,
            {"epsilon": 1.2, "sigma": 2.2, "cutoff": 3.2, "shift": 4.2,
             "offset": 5.2, "e1": 72, "e2": 82, "b1": 9.2, "b2": 10.2},
            "generic_lennard_jones")

    if espressomd.has_features(["GAY_BERNE"]):
        test_gb = generateTestForNon_bonded_interaction(
            0, 0, espressomd.interactions.GayBerneInteraction,
            {"eps": 1.0, "sig": 1.0, "cut": 4.0, "k1": 3.0,
             "k2": 5.0, "mu": 2.0, "nu": 1.0},
            "gay_berne")

    @utx.skipIfMissingFeatures(["LENNARD_JONES", "WCA"])
    def test_set_params(self):
        self.system.non_bonded_inter[0, 0].lennard_jones.set_params(
            epsilon=1., sigma=2., cutoff=3., shift="auto")
        self.system.non_bonded_inter.reset()
        self.assertEqual(self.system.non_bonded_inter[0, 0].lennard_jones.get_params(),
                         {'shift': 0., 'sigma': 0., 'epsilon': 0.,
                          'cutoff': -1., 'offset': 0., 'min': 0.})
        wca = espressomd.interactions.WCAInteraction(epsilon=1., sigma=2.)
        wca.set_params(epsilon=2.)
        self.assertEqual(wca.get_params()["epsilon"], 2.)
        self.system.non_bonded_inter[0, 0].wca = wca
        self.assertEqual(
            self.system.non_bonded_inter[0, 0].wca.get_params()["epsilon"], 2.)
        self.assertEqual(wca.get_params()["epsilon"], 2.)
        wca.set_params(epsilon=3.)
        self.assertEqual(
            self.system.non_bonded_inter[0, 0].wca.get_params()["epsilon"], 3.)
        self.assertEqual(wca.get_params()["epsilon"], 3.)
        self.system.non_bonded_inter.reset()
        wca.set_params(epsilon=4.)
        self.assertEqual(
            self.system.non_bonded_inter[0, 0].wca.get_params()["epsilon"], 4.)
        self.assertEqual(wca.get_params()["epsilon"], 4.)
        with self.assertRaisesRegex(RuntimeError, r"Non-bonded interaction is already bound to interaction pair \[0, 0\]"):
            self.system.non_bonded_inter[0, 1].wca = wca

    @utx.skipIfMissingFeatures("LENNARD_JONES")
    def test_exceptions(self):
        err_msg_required = (r"The following keys have to be given as keyword arguments: "
                            r"\['cutoff', 'epsilon', 'shift', 'sigma'\], got "
                            r"\['epsilon', 'sigma'\] \(missing \['cutoff', 'shift'\]\)")
        err_msg_valid = (r"Only the following keys can be given as keyword arguments: "
                         r"\['cutoff', 'epsilon', 'min', 'offset', 'shift', 'sigma'\], got "
                         r"\['cutoff', 'epsilon', 'shift', 'sigma', 'unknown'\] \(unknown \['unknown'\]\)")
        with self.assertRaisesRegex(ValueError, err_msg_required):
            espressomd.interactions.LennardJonesInteraction(
                epsilon=1., sigma=2.)
        with self.assertRaisesRegex(ValueError, err_msg_required):
            self.system.non_bonded_inter[0, 0].lennard_jones.set_params(
                epsilon=1., sigma=2.)
        with self.assertRaisesRegex(ValueError, err_msg_valid):
            espressomd.interactions.LennardJonesInteraction(
                epsilon=1., sigma=2., cutoff=3., shift=4., unknown=5.)
        with self.assertRaisesRegex(ValueError, err_msg_valid):
            self.system.non_bonded_inter[0, 0].lennard_jones.set_params(
                epsilon=1., sigma=2., cutoff=3., shift=4., unknown=5.)
        with self.assertRaisesRegex(ValueError, "LJ parameter 'shift' has to be 'auto' or a float"):
            self.system.non_bonded_inter[0, 0].lennard_jones.set_params(
                epsilon=1., sigma=2., cutoff=3., shift="automatic")

        skin = self.system.cell_system.skin
        box_l = self.system.box_l
        node_grid = self.system.cell_system.node_grid
        n_nodes = self.system.cell_system.get_state()['n_nodes']
        max_ia_cutoff = min(box_l / node_grid) - skin * (n_nodes > 1)
        wrong_cutoff = 1.01 * max_ia_cutoff
        lennard_jones = self.system.non_bonded_inter[0, 0].lennard_jones
        with self.assertRaisesRegex(Exception, "setting LennardJones raised an error: ERROR: interaction range .+ in direction [0-2] is larger than the local box size"):
            lennard_jones.set_params(
                epsilon=1., sigma=1., cutoff=wrong_cutoff, shift="auto")
        self.assertAlmostEqual(
            lennard_jones.get_params()['cutoff'], wrong_cutoff, delta=1e-10)

    def check_potential_exceptions(self, ia_class, params, check_keys):
        for key in check_keys:
            with self.assertRaisesRegex(ValueError, f"parameter '{key}'"):
                invalid_params = params.copy()
                invalid_params[key] = -0.1
                ia_class(**invalid_params)

    @utx.skipIfMissingFeatures("WCA")
    def test_wca_exceptions(self):
        self.check_potential_exceptions(
            espressomd.interactions.WCAInteraction,
            {"epsilon": 1., "sigma": 1.},
            ("epsilon", "sigma"))

    @utx.skipIfMissingFeatures("LENNARD_JONES")
    def test_lj_exceptions(self):
        self.check_potential_exceptions(
            espressomd.interactions.LennardJonesInteraction,
            {"epsilon": 1., "sigma": 1., "cutoff": 1.5, "shift": 0.2},
            ("epsilon", "sigma"))

    @utx.skipIfMissingFeatures("LJCOS")
    def test_lj_cos_exceptions(self):
        self.check_potential_exceptions(
            espressomd.interactions.LennardJonesCosInteraction,
            {"epsilon": 1., "sigma": 1., "cutoff": 1.5, "offset": 0.2},
            ("epsilon", "sigma"))

    @utx.skipIfMissingFeatures("LJCOS2")
    def test_lj_cos2_exceptions(self):
        self.check_potential_exceptions(
            espressomd.interactions.LennardJonesCos2Interaction,
            {"epsilon": 1., "sigma": 1., "width": 1.5, "offset": 0.2},
            ("epsilon", "sigma"))

    @utx.skipIfMissingFeatures("HERTZIAN")
    def test_hertzian_exceptions(self):
        self.check_potential_exceptions(
            espressomd.interactions.HertzianInteraction,
            {"eps": 1., "sig": 1.},
            ("eps",)
        )

    @utx.skipIfMissingFeatures("GAUSSIAN")
    def test_gaussian_exceptions(self):
        self.check_potential_exceptions(
            espressomd.interactions.GaussianInteraction,
            {"eps": 1., "sig": 1., "cutoff": 1.5},
            ("eps", "sig")
        )


if __name__ == "__main__":
    ut.main()
