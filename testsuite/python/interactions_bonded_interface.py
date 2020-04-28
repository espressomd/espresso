#
# Copyright (C) 2013-2019 The ESPResSo project
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
from tests_common import assert_params_match

import espressomd


class ParticleProperties(ut.TestCase):
    system = espressomd.System(box_l=[20.0, 20.0, 20.0])

    # Particle id to work on
    pid = 17

    def bondsMatch(self, inType, outBond, inParams, outParams, msg_long):
        """Check, if the bond type set and gotten back as well as the bond
        parameters set and gotten back match. Only check keys present in
        ``inParams``.
        """
        self.assertEqual(outBond, inType, msg="Bonded interaction mismatch")
        assert_params_match(self, inParams, outParams, msg_long)

    def parameterKeys(self, bondObject):
        """
        Check :meth:`~espressomd.interactions.NonBondedInteraction.valid_keys`
        and :meth:`~espressomd.interactions.NonBondedInteraction.required_keys`
        return sets, and that
        :meth:`~espressomd.interactions.NonBondedInteraction.set_default_params`
        returns a dictionary with the correct keys.

        Parameters
        ----------
        bondObject: instance of a class derived from :class:`espressomd.interactions.BondedInteraction`
            Object of the interaction to test, e.g.
            :class:`~espressomd.interactions.FeneBond`
        """
        classname = bondObject.__class__.__name__
        valid_keys = bondObject.valid_keys()
        required_keys = bondObject.required_keys()
        old_params = dict(bondObject.params)
        bondObject.set_default_params()
        default_keys = set(bondObject.params.keys())
        bondObject.params = old_params
        self.assertIsInstance(valid_keys, set,
                              "{}.valid_keys() must return a set".format(
                                  classname))
        self.assertIsInstance(required_keys, set,
                              "{}.required_keys() must return a set".format(
                                  classname))
        self.assertTrue(default_keys.issubset(valid_keys),
                        "{}.set_default_params() has unknown parameters: {}".format(
            classname, default_keys.difference(valid_keys)))
        self.assertTrue(default_keys.isdisjoint(required_keys),
                        "{}.set_default_params() has extra parameters: {}".format(
            classname, default_keys.intersection(required_keys)))
        self.assertSetEqual(default_keys, valid_keys - required_keys,
                            "{}.set_default_params() should have keys: {}, got: {}".format(
                                classname, valid_keys - required_keys, default_keys))

    def setUp(self):
        if not self.system.part.exists(self.pid):
            self.system.part.add(id=self.pid, pos=(0, 0, 0, 0))

    def generateTestForBondParams(_bondId, _bondClass, _params):
        """Generates test cases for checking bond parameters set and gotten
        back from the espresso core actually match those in the Python classes.
        Only keys which are present in ``_params`` are checked.

        Parameters
        ----------
        _bondId: :obj:`int`
            Identifier of the bonded IA in Espresso to test on, e.g. 0, 2, 1...
        _bondClass: class derived from :class:`espressomd.interactions.BondedInteraction`
            Class of the interaction to test, e.g.
            :class:`~espressomd.interactions.FeneBond`
        _params: :obj:`dict`
            Bond parameters, e.g. ``{"k": 1., "r_0": 0}``
        """
        bondId = _bondId
        bondClass = _bondClass
        params = _params

        def func(self):
            # This code is run at the execution of the generated function.
            # It will use the state of the variables in the outer function,
            # which was there, when the outer function was called
            self.system.bonded_inter[bondId] = bondClass(**params)
            outBond = self.system.bonded_inter[bondId]
            tnIn = bondClass(**params).type_number()
            tnOut = outBond.type_number()
            outParams = outBond.params
            self.bondsMatch(
                tnIn, tnOut, params, outParams,
                "{}: value set and value gotten back differ for bond id {}: {} vs. {}"
                .format(bondClass(**params).type_name(), bondId, params, outParams))
            self.parameterKeys(outBond)

        return func

    test_fene = generateTestForBondParams(
        0, espressomd.interactions.FeneBond, {"r_0": 1.1, "k": 5.2, "d_r_max": 3.})
    test_fene2 = generateTestForBondParams(
        1, espressomd.interactions.FeneBond, {"r_0": 1.1, "k": 5.2, "d_r_max": 3.})
    test_harmonic = generateTestForBondParams(
        0, espressomd.interactions.HarmonicBond, {"r_0": 1.1, "k": 5.2})
    test_harmonic2 = generateTestForBondParams(
        0, espressomd.interactions.HarmonicBond, {"r_0": 1.1, "k": 5.2, "r_cut": 1.3})

    # HarmonicDumbbell has only interface tests, so it is marked as
    # experimental
    if espressomd.has_features(["ROTATION", "EXPERIMENTAL_FEATURES"]):
        test_harmonic_dumbbell = generateTestForBondParams(
            0, espressomd.interactions.HarmonicDumbbellBond, {"k1": 1.1, "k2": 2.2, "r_0": 1.5})
        test_harmonic_dumbbell2 = generateTestForBondParams(
            0, espressomd.interactions.HarmonicDumbbellBond, {"k1": 1.1, "k2": 2.2, "r_0": 1.5, "r_cut": 1.9})

    test_dihedral = generateTestForBondParams(
        0, espressomd.interactions.Dihedral, {"mult": 3.0, "bend": 5.2, "phase": 3.})

    test_angle_harm = generateTestForBondParams(
        0, espressomd.interactions.AngleHarmonic, {"bend": 5.2, "phi0": 3.2})
    test_angle_cos = generateTestForBondParams(
        0, espressomd.interactions.AngleCosine, {"bend": 5.2, "phi0": 3.2})
    test_angle_cossquare = generateTestForBondParams(
        0, espressomd.interactions.AngleCossquare, {"bend": 5.2, "phi0": 0.})

    test_tabulated_bond = generateTestForBondParams(
        0, espressomd.interactions.TabulatedDistance, {"min": 1.,
                                                       "max": 2.,
                                                       "energy": [1., 2., 3.],
                                                       "force": [3., 4., 5.]})
    test_tabulated = generateTestForBondParams(
        0, espressomd.interactions.TabulatedAngle, {"energy": [1., 2., 3.],
                                                    "force": [3., 4., 5.]})
    test_tabulated = generateTestForBondParams(
        0, espressomd.interactions.TabulatedDihedral, {"energy": [1., 2., 3.],
                                                       "force": [3., 4., 5.]})


if __name__ == "__main__":
    ut.main()
