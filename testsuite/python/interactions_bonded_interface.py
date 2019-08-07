#
# Copyright (C) 2013-2018 The ESPResSo project
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
import espressomd


class ParticleProperties(ut.TestCase):
    system = espressomd.System(box_l=[10.0, 10.0, 10.0])

    # Particle id to work on
    pid = 17

    # Error tolerance when comparing arrays/tuples...
    tol = 1E-9

    def bondsMatch(self, inType, outType, inParams, outParams):
        """Check, if the bond type set and gotten back as well as the bond
        parameters set and gotten back match. Only check keys present in
        inParams.
        """
        if inType != outType:
            return False

        for k in list(inParams.keys()):
            if k not in outParams:
                return False
            if outParams[k] != inParams[k]:
                return False

        return True

    def setUp(self):
        if not self.system.part.exists(self.pid):
            self.system.part.add(id=self.pid, pos=(0, 0, 0, 0))

    def generateTestForBondParams(_bondId, _bondClass, _params):
        """Generates test cases for checking bond parameters set and gotten back
        from Es actually match. Only keys which are present in _params are checked
        1st arg: Id of the bonded ia in Espresso to test on, i.e., 0,2,1...
        2nd: Class of the bond potential to test, ie.e, FeneBond, HarmonicBond
        3rd: Bond parameters as dictionary, i.e., {"k"=1.,"r_0"=0}
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
            self.assertTrue(
                self.bondsMatch(
                    tnIn,
                    tnOut,
                    params,
                    outParams),
                bondClass(**params).type_name() +
                ": value set and value gotten back differ for bond id " +
                str(bondId) +
                ": " +
                params.__str__() +
                " vs. " +
                outParams.__str__())

        return func

    test_fene = generateTestForBondParams(
        0, espressomd.interactions.FeneBond, {"r_0": 1.1, "k": 5.2, "d_r_max": 3.})
    test_fene2 = generateTestForBondParams(
        1, espressomd.interactions.FeneBond, {"r_0": 1.1, "k": 5.2, "d_r_max": 3.})
    test_harmonic = generateTestForBondParams(
        0, espressomd.interactions.HarmonicBond, {"r_0": 1.1, "k": 5.2})
    test_harmonic2 = generateTestForBondParams(
        0, espressomd.interactions.HarmonicBond, {"r_0": 1.1, "k": 5.2, "r_cut": 1.3})

    if espressomd.has_features(["ROTATION"]):
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
    if espressomd.has_features(["LENNARD_JONES"]):
        test_subt_lj = generateTestForBondParams(
            0, espressomd.interactions.SubtLJ, {})

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
