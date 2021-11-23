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
import unittest_decorators as utx
import tests_common

import espressomd
import espressomd.interactions


class BondedInteractions(ut.TestCase):
    system = espressomd.System(box_l=[20.0, 20.0, 20.0])

    def setUp(self):
        self.system.part.add(pos=4 * [[0, 0, 0]])

    def tearDown(self):
        self.system.part.clear()
        self.system.bonded_inter.clear()

    def bondsMatch(self, inType, outBond, inParams, outParams, msg_long):
        """Check, if the bond type set and gotten back as well as the bond
        parameters set and gotten back match. Only check keys present in
        ``inParams``.
        """
        self.assertEqual(outBond, inType, msg="Bonded interaction mismatch")
        tests_common.assert_params_match(self, inParams, outParams, msg_long)

    def parameterKeys(self, bondObject):
        """
        Check :meth:`~espressomd.interactions.BondedInteraction.valid_keys`
        and :meth:`~espressomd.interactions.BondedInteraction.required_keys`
        return sets, and that
        :meth:`~espressomd.interactions.BondedInteraction.get_default_params`
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
        old_params = bondObject.params.copy()
        default_keys = set(bondObject.get_default_params())
        self.assertIsInstance(valid_keys, set,
                              f"{classname}.valid_keys() must return a set")
        self.assertIsInstance(required_keys, set,
                              f"{classname}.required_keys() must return a set")
        self.assertTrue(default_keys.issubset(valid_keys),
                        f"{classname}.get_default_params() has unknown "
                        f"parameters: {default_keys.difference(valid_keys)}")
        self.assertTrue(default_keys.isdisjoint(required_keys),
                        f"{classname}.get_default_params() has extra "
                        f"parameters: {default_keys.intersection(required_keys)}")
        self.assertSetEqual(default_keys, valid_keys - required_keys,
                            f"{classname}.get_default_params() should have keys: "
                            f"{valid_keys - required_keys}, got: {default_keys}")
        with self.assertRaisesRegex(RuntimeError, "Bond parameters are immutable"):
            bondObject.params = {}
        self.assertEqual(bondObject.params, old_params)

    def generateTestForBondParams(_bondId, _bondClass, _params, _refs=None):
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
        _refs: :obj:`set` or :obj:`dict`
            Subset of keys in ``_params`` to check in the bond parameters read
            from the core, or subset of the ``_params`` dictionary. Leave it
            as ``None`` if all keys need to be checked.
        """
        bondId = _bondId
        bondClass = _bondClass
        params = _params
        if _refs is None:
            outParamsRef = _params.copy()
        elif isinstance(_refs, dict):
            outParamsRef = _refs
        else:
            outParamsRef = {key: params[key] for key in _refs}

        def func(self):
            # This code is run at the execution of the generated function.
            # It will use the state of the variables in the outer function,
            # which was there, when the outer function was called
            self.system.bonded_inter.remove(bondId)
            self.system.bonded_inter.remove(bondId + 1)

            bond = bondClass(**params)
            self.system.bonded_inter[bondId] = bond
            self.system.bonded_inter[bondId + 1] = bond

            # check that bonds are identified by their pointer
            self.assertEqual(bond, self.system.bonded_inter[bondId])
            self.assertEqual(bond, self.system.bonded_inter[bondId + 1])
            # check that bonds are removed
            bonded_len = len(self.system.bonded_inter)
            self.system.bonded_inter.remove(bondId + 1)
            new_bonded_len = len(self.system.bonded_inter)
            self.assertEqual(new_bonded_len, bonded_len - 1)
            with self.assertRaises(Exception):
                self.system.bonded_inter[bondId + 1]
            # check that bonds are distinguished by their internal object address
            # put two identical bonds in the same ID, so that they only differ
            # in their shared_ptr address
            new_bond = bondClass(**params)
            self.system.bonded_inter[bondId] = new_bond
            self.system.bonded_inter[bondId] = bond
            self.assertNotEqual(bond, new_bond)

            # check that parameters written and read back are identical
            outBond = self.system.bonded_inter[bondId]
            tnIn = bondClass(**params).type_number()
            tnOut = outBond.type_number()
            outParams = outBond.params
            self.bondsMatch(
                tnIn, tnOut, outParamsRef, outParams,
                "{}: value set and value gotten back differ for bond id {}: {} vs. {}"
                .format(bondClass(**params).type_name(), bondId, outParamsRef, outParams))
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
    test_quartic = generateTestForBondParams(
        0, espressomd.interactions.QuarticBond, {"k0": 2., "k1": 5., "r": 0.5, "r_cut": 1.2})
    test_ibm_volcons = generateTestForBondParams(
        0, espressomd.interactions.IBM_VolCons, {"softID": 15, "kappaV": 0.01})
    test_ibm_tribend = generateTestForBondParams(
        0, espressomd.interactions.IBM_Tribend,
        {"ind1": 0, "ind2": 1, "ind3": 2, "ind4": 3,
            "kb": 1.1, "refShape": "Initial"},
        {"kb": 1.1, "theta0": 0.0})
    test_ibm_triel = generateTestForBondParams(
        0, espressomd.interactions.IBM_Triel,
        {"ind1": 0, "ind2": 1, "ind3": 2, "k1": 1.1, "k2": 1.2,
            "maxDist": 1.6, "elasticLaw": "NeoHookean"},
        {"k1", "k2", "maxDist", "elasticLaw"})

    test_dihedral = generateTestForBondParams(
        0, espressomd.interactions.Dihedral, {"mult": 3, "bend": 5.2, "phase": 3.})

    test_angle_harm = generateTestForBondParams(
        0, espressomd.interactions.AngleHarmonic, {"bend": 5.2, "phi0": 3.2})
    test_angle_cos = generateTestForBondParams(
        0, espressomd.interactions.AngleCosine, {"bend": 5.2, "phi0": 3.2})
    test_angle_cossquare = generateTestForBondParams(
        0, espressomd.interactions.AngleCossquare, {"bend": 5.2, "phi0": 0.})

    @utx.skipIfMissingFeatures(["TABULATED"])
    def test_tabulated_bond(self):
        params = {"min": 1., "max": 2., "energy": [1., 2., 3.],
                  "force": [3., 4., 5.]}
        BondedInteractions.generateTestForBondParams(
            0, espressomd.interactions.TabulatedDistance, params)(self)

    @utx.skipIfMissingFeatures(["TABULATED"])
    def test_tabulated_angle(self):
        params = {"energy": [1., 2., 3.], "force": [3., 4., 5.]}
        BondedInteractions.generateTestForBondParams(
            0, espressomd.interactions.TabulatedAngle, params)(self)

    @utx.skipIfMissingFeatures(["TABULATED"])
    def test_tabulated_dihedral(self):
        params = {"energy": [1., 2., 3.], "force": [3., 4., 5.]}
        BondedInteractions.generateTestForBondParams(
            0, espressomd.interactions.TabulatedDihedral, params)(self)

    @utx.skipIfMissingFeatures(["BOND_CONSTRAINT"])
    def test_rigid_bond(self):
        params = {"r": 1.2, "ptol": 1E-3, "vtol": 1E-3}
        BondedInteractions.generateTestForBondParams(
            2, espressomd.interactions.RigidBond, params)(self)

    @utx.skipIfMissingFeatures(["ELECTROSTATICS"])
    def test_bonded_coulomb(self):
        params = {"prefactor": 2.46}
        BondedInteractions.generateTestForBondParams(
            0, espressomd.interactions.BondedCoulomb, params)(self)

    @utx.skipIfMissingFeatures(["ELECTROSTATICS"])
    def test_bonded_coulomb_sr(self):
        params = {"q1q2": 3.46}
        BondedInteractions.generateTestForBondParams(
            0, espressomd.interactions.BondedCoulombSRBond, params)(self)

    def test_exceptions(self):
        error_msg_not_yet_defined = 'The bond with id 0 is not yet defined'
        bond_type = espressomd.interactions.get_bonded_interaction_type_from_es_core(
            5000)
        self.assertEqual(bond_type, 0)
        has_bond = self.system.bonded_inter.call_method('has_bond', bond_id=0)
        self.assertFalse(has_bond)
        with self.assertRaisesRegex(ValueError, error_msg_not_yet_defined):
            self.system.bonded_inter[0]
        with self.assertRaisesRegex(IndexError, error_msg_not_yet_defined):
            self.system.bonded_inter.call_method('get_bond', bond_id=0)

        # bonds can only be overwritten by bonds of the same type
        harm_bond1 = espressomd.interactions.HarmonicBond(r_0=1., k=1.)
        harm_bond2 = espressomd.interactions.HarmonicBond(r_0=2., k=2.)
        fene_bond = espressomd.interactions.FeneBond(r_0=1., k=1., d_r_max=2.)
        angle_bond = espressomd.interactions.AngleHarmonic(bend=5., phi0=3.)
        self.system.bonded_inter[0] = harm_bond1
        self.system.bonded_inter[0] = harm_bond2
        with self.assertRaisesRegex(ValueError, 'Bonds can only be overwritten by bonds of equal type'):
            self.system.bonded_inter[0] = fene_bond
        with self.assertRaisesRegex(ValueError, 'Bonds can only be overwritten by bonds of equal type'):
            self.system.bonded_inter[0] = angle_bond

        # bonds are immutable
        with self.assertRaisesRegex(RuntimeError, "Parameter 'r_0' is read-only"):
            harm_bond1.r_0 = 5.

        # sanity checks during bond construction
        with self.assertRaisesRegex(RuntimeError, "Parameter 'r_0' is missing"):
            espressomd.interactions.HarmonicBond(k=1.)
        with self.assertRaisesRegex(ValueError, "Unknown refShape: 'Unknown'"):
            espressomd.interactions.IBM_Tribend(
                ind1=0, ind2=1, ind3=2, ind4=3, kb=1.1, refShape='Unknown')
        with self.assertRaisesRegex(ValueError, "Unknown elasticLaw: 'Unknown'"):
            espressomd.interactions.IBM_Triel(
                ind1=0, ind2=1, ind3=2, k1=1.1, k2=1.2, maxDist=1.6, elasticLaw='Unknown')

        # sanity checks when removing bonds
        self.system.bonded_inter.clear()
        with self.assertRaisesRegex(ValueError, error_msg_not_yet_defined):
            self.system.bonded_inter[0]
        self.system.bonded_inter[0] = harm_bond1
        self.system.bonded_inter[0]
        self.system.bonded_inter.remove(0)
        with self.assertRaisesRegex(ValueError, error_msg_not_yet_defined):
            self.system.bonded_inter[0]
        self.system.bonded_inter[0] = harm_bond1
        self.system.bonded_inter[0]
        del self.system.bonded_inter[0]
        with self.assertRaisesRegex(ValueError, error_msg_not_yet_defined):
            self.system.bonded_inter[0]


if __name__ == "__main__":
    ut.main()
