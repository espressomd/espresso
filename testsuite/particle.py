#
# Copyright (C) 2013,2014,2015,2016 The ESPResSo project
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
# Tests particle property setters/getters
from __future__ import print_function
import unittest as ut
import espressomd
import numpy as np
from espressomd.interactions import FeneBond


class ParticleProperties(ut.TestCase):

    # Particle id to work on
    pid = 17

    # Error tolerance when comparing arrays/tuples...
    tol = 1E-9

    # Handle for espresso system
    es = espressomd.System()

    f1 = FeneBond(k=1, d_r_max=5)
    es.bonded_inter.add(f1)
    f2 = FeneBond(k=1, d_r_max=5)
    es.bonded_inter.add(f2)

    def arraysNearlyEqual(self, a, b):
        """Test, if the magnitude of the difference between two arrays is smaller than the tolerance"""

        # Check length
        if len(a) != len(b):
            return False

        # We have to use a loop, since we can't be sure, we're getting numpy
        # arrays
        sum = 0.
        for i in range(len(a)):
            sum += abs(a[i] - b[i])

        if sum > self.tol:
            return False

        return True

    def setUp(self):
        if not self.es.part.exists(self.pid):
            self.es.part.add(id=self.pid, pos=(0, 0, 0))

    def generateTestForVectorProperty(_propName, _value):
        """Generates test cases for vectorial particle properties such as
        position, velocity...
        1st arg: name of the property (e.g., "pos"),
        2nd array: value to be used for testing. Has to be numpy.array of floats
        """
        # This is executed, when generateTestForVectorProperty() is called
        propName = _propName
        value = _value

        def func(self):
            # This code is run at the execution of the generated function.
            # It will use the state of the variables in the outer function,
            # which was there, when the outer function was called
            setattr(self.es.part[self.pid], propName, value)
            self.assertTrue(self.arraysNearlyEqual(getattr(self.es.part[
                            self.pid], propName), value), propName + ": value set and value gotten back differ.")

        return func

    def generateTestForScalarProperty(_propName, _value):
        """Generates test cases for scalar particle properties such as
        type, mass, charge...
        1st arg: name of the property (e.g., "type"),
        2nd array: value to be used for testing. int or float
        """
        # This is executed, when generateTestForVectorProperty() is called
        propName = _propName
        value = _value

        def func(self):
            # This code is run at the execution of the generated function.
            # It will use the state of the variables in the outer function,
            # which was there, when the outer function was called
            setattr(self.es.part[self.pid], propName, value)
            self.assertEqual(getattr(self.es.part[self.pid], propName),
                             value, propName + ": value set and value gotten back differ.")

        return func

    test_pos = generateTestForVectorProperty("pos", np.array([0.1, 0.2, 0.3]))
    test_v = generateTestForVectorProperty("v", np.array([0.2, 0.3, 0.4]))
    test_f = generateTestForVectorProperty("f", np.array([0.2, 0.3, 0.7]))
    test_type = generateTestForScalarProperty("type", int(3))
    test_mol_id = generateTestForScalarProperty("mol_id", int(3))

    test_bonds_property = generateTestForScalarProperty(
        "bonds", ((f1, 1), (f2, 2)))

    if espressomd.has_features(["MASS"]):
        test_mass = generateTestForScalarProperty("mass", 1.3)

    if espressomd.has_features(["ROTATION"]):
        
        for x in 0,1:
            for y in 0,1:
                for z in 0,1:
                    test_rotation = generateTestForVectorProperty(
                        "rotation", np.array([x,y,z],dtype=int))

        test_omega_lab = generateTestForVectorProperty(
            "omega_lab", np.array([4., 2., 1.]))
        test_omega_body = generateTestForVectorProperty(
            "omega_body", np.array([4., 72., 1.]))
        test_torque_lab = generateTestForVectorProperty(
            "torque_lab", np.array([4., 72., 3.7]))
        # The tested value has to be normalized!
        test_quat = generateTestForVectorProperty(
            "quat", np.array([0.5, 0.5, 0.5, 0.5]))

        if espressomd.has_features(["LANGEVIN_PER_PARTICLE"]):
            if espressomd.has_features(["PARTICLE_ANISOTROPY"]):
                test_gamma = generateTestForVectorProperty(
                    "gamma", np.array([2., 9., 0.23]))
            else:
                test_gamma = generateTestForScalarProperty("gamma", 17.3)

            if espressomd.has_features(["ROTATIONAL_INERTIA"]):
                test_gamma_rot = generateTestForVectorProperty(
                    "gamma_rot", np.array([5., 10., 0.33]))
            else:
                test_gamma_rot = generateTestForScalarProperty(
                    "gamma_rot", 14.23)
#    test_director=generateTestForVectorProperty("director",np.array([0.5,0.4,0.3]))

    if espressomd.has_features(["ELECTROSTATICS"]):
        test_charge = generateTestForScalarProperty("q", -19.7)

    if espressomd.has_features(["DIPOLES"]):
        test_dip = generateTestForVectorProperty(
            "dip", np.array([0.5, -0.5, 3]))
        test_dipm = generateTestForScalarProperty("dipm", -9.7)

    if espressomd.has_features(["VIRTUAL_SITES"]):
        test_virtual = generateTestForScalarProperty("virtual", 1)
    if espressomd.has_features(["VIRTUAL_SITES_RELATIVE"]):
        def test_zz_vs_relative(self):
            self.es.part.add(id=0, pos=(0, 0, 0))
            self.es.part.add(id=1, pos=(0, 0, 0))
            self.es.part[1].vs_relative = (0, 5.0, (0.5, -0.5, -0.5, -0.5))
            res = self.es.part[1].vs_relative
            self.assertEqual(res[0], 0, "vs_relative: " + res.__str__())
            self.assertEqual(res[1], 5.0, "vs_relative: " + res.__str__())
            self.assertTrue(self.arraysNearlyEqual(res[2], np.array(
                (0.5, -0.5, -0.5, -0.5))), "vs_relative: " + res.__str__())


if __name__ == "__main__":
    #print("Features: ", espressomd.features())
    ut.main()
