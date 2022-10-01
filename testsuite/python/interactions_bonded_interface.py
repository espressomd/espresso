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
import espressomd.code_features
import numpy as np


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
            tnIn = bondClass._type_number
            tnOut = outBond._type_number
            outParams = outBond.params
            self.bondsMatch(
                tnIn, tnOut, outParamsRef, outParams,
                f"{bondClass.__name__}: value set and value gotten back "
                f"differ for bond id {bondId}: {outParamsRef} vs. {outParams}")
            with self.assertRaisesRegex(RuntimeError, "Bond parameters are immutable"):
                outBond.params = {}
            old_params = outBond.params.copy()
            for key in (outBond.params.keys() | old_params.keys()):
                if isinstance(old_params[key], str):
                    self.assertEqual(outBond.params[key], old_params[key])
                else:
                    np.testing.assert_allclose(
                        outBond.params[key], old_params[key], atol=1e-10)

            # check no-op
            self.assertIsNone(outBond.call_method('unknown'))

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
    test_ibm_tribend_flat = generateTestForBondParams(
        0, espressomd.interactions.IBM_Tribend,
        {"ind1": 0, "ind2": 1, "ind3": 2, "ind4": 3,
            "kb": 1.1, "refShape": "Flat"},
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
        bond_type = self.system.bonded_inter.call_method(
            'get_zero_based_type', bond_id=5000)
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
        with self.assertRaisesRegex(RuntimeError, "No bond with id 8 exists in the ESPResSo core"):
            espressomd.interactions.FeneBond(bond_id=8)
        with self.assertRaisesRegex(RuntimeError, "The bond with id 0 is not defined as a FENE bond in the ESPResSo core"):
            espressomd.interactions.FeneBond(bond_id=0)

        # bonds can only be compared for equality
        self.assertEqual(angle_bond, angle_bond)
        self.assertNotEqual(harm_bond1, harm_bond2)
        self.assertNotEqual(angle_bond, fene_bond)
        with self.assertRaises(NotImplementedError):
            angle_bond > fene_bond
        with self.assertRaises(NotImplementedError):
            angle_bond <= fene_bond

        # bonds are immutable
        with self.assertRaisesRegex(RuntimeError, "Parameter 'r_0' is read-only"):
            harm_bond1.r_0 = 5.

        # sanity checks during bond construction
        with self.assertRaisesRegex(RuntimeError, "Parameter 'r_0' is missing"):
            espressomd.interactions.HarmonicBond(k=1.)
        with self.assertRaisesRegex(RuntimeError, f"Parameter 'rcut' is not recognized"):
            espressomd.interactions.HarmonicBond(k=1., r_0=1., rcut=2.)
        with self.assertRaisesRegex(ValueError, "Invalid value for parameter 'refShape': 'Unknown'"):
            espressomd.interactions.IBM_Tribend(
                ind1=0, ind2=1, ind3=2, ind4=3, kb=1.1, refShape='Unknown')
        with self.assertRaisesRegex(ValueError, "Invalid value for parameter 'elasticLaw': 'Unknown'"):
            espressomd.interactions.IBM_Triel(
                ind1=0, ind2=1, ind3=2, k1=1.1, k2=1.2, maxDist=1.6, elasticLaw='Unknown')
        with self.assertRaisesRegex(ValueError, "A parameter 'seed' has to be given on first activation of a thermalized bond"):
            espressomd.interactions.ThermalizedBond(
                temp_com=1., gamma_com=1., temp_distance=1., gamma_distance=1.,
                r_cut=2.)
        with self.assertRaisesRegex(ValueError, "Parameter 'seed' must be >= 0"):
            espressomd.interactions.ThermalizedBond(
                temp_com=1., gamma_com=1., temp_distance=1., gamma_distance=1.,
                r_cut=2., seed=-1)

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

        # check cutoff exceptions
        skin = self.system.cell_system.skin
        box_l = self.system.box_l
        node_grid = self.system.cell_system.node_grid
        n_nodes = self.system.cell_system.get_state()['n_nodes']
        max_ia_cutoff = min(box_l / node_grid) - skin * (n_nodes > 1)
        err_msg = r"while calling method insert\(\): ERROR: "
        if n_nodes == 1:
            safe_cut = 0.5 * max_ia_cutoff
            err_msg += r"number of cells 1 is smaller than minimum 8"
        else:
            safe_cut = 1.0 * max_ia_cutoff
            err_msg += r"interaction range .+ in direction [0-2] is larger than the local box size"

        # a bond with exactly the right cutoff can be added
        h0 = espressomd.interactions.HarmonicBond(k=1., r_0=0., r_cut=safe_cut)
        self.system.bonded_inter.add(h0)
        self.assertEqual(len(self.system.bonded_inter), 1)
        self.system.bonded_inter.clear()

        # On more than one node:
        # a bond with a cutoff larger than the box can be added, but throws
        big_cut = 1.001 * safe_cut
        if n_nodes > 1:
            h1 = espressomd.interactions.HarmonicBond(
                k=1., r_0=0., r_cut=big_cut)
            with self.assertRaisesRegex(Exception, err_msg):
                self.system.bonded_inter.add(h1)
            self.assertEqual(len(self.system.bonded_inter), 1)
            self.system.bonded_inter.clear()

        # a dihedral halves the cutoff
        safe_cut /= 2.
        h2 = espressomd.interactions.HarmonicBond(k=1., r_0=0., r_cut=safe_cut)
        self.system.bonded_inter.add(h2)
        dihe = espressomd.interactions.Dihedral(bend=1., mult=1, phase=0.)
        self.system.bonded_inter.add(dihe)
        self.assertEqual(len(self.system.bonded_inter), 2)
        self.system.bonded_inter.clear()

        # a dihedral halves the cutoff, safe cutoffs become unsafe
        if n_nodes > 1:
            half_cut = big_cut / 2.
            h3 = espressomd.interactions.HarmonicBond(
                k=1., r_0=0., r_cut=half_cut)
            self.system.bonded_inter.add(h3)
            with self.assertRaisesRegex(Exception, err_msg):
                self.system.bonded_inter.add(dihe)
            self.assertEqual(len(self.system.bonded_inter), 2)
            self.system.bonded_inter.clear()

    def test_feature_checks(self):
        base_class = espressomd.interactions.BondedInteraction
        FeaturesError = espressomd.code_features.FeaturesError
        for class_bond in espressomd.interactions.__dict__.values():
            if isinstance(class_bond, type) and issubclass(
                    class_bond, base_class) and class_bond != base_class:
                feature = getattr(class_bond, "_so_feature", None)
                if feature is not None and not espressomd.code_features.has_features(
                        feature):
                    with self.assertRaisesRegex(FeaturesError, f"Missing features {feature}"):
                        class_bond()


if __name__ == "__main__":
    ut.main()
