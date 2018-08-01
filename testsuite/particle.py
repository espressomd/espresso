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
    system = espressomd.System(box_l=[100.0, 100.0, 100.0])

    f1 = FeneBond(k=1, d_r_max=5)
    system.bonded_inter.add(f1)
    f2 = FeneBond(k=1, d_r_max=5)
    system.bonded_inter.add(f2)

    def setUp(self):
        if not self.system.part.exists(self.pid):
            self.system.part.add(id=self.pid, pos=(0, 0, 0))

    def tearDown(self):
        self.system.part.clear()

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
            setattr(self.system.part[self.pid], propName, value)
            np.testing.assert_allclose(np.array(getattr(self.system.part[
                self.pid], propName)), value, err_msg=propName + ": value set and value gotten back differ.", atol=self.tol)

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
            setattr(self.system.part[self.pid], propName, value)
            self.assertEqual(getattr(self.system.part[self.pid], propName),
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

        for x in 0, 1:
            for y in 0, 1:
                for z in 0, 1:
                    test_rotation = generateTestForVectorProperty(
                        "rotation", np.array([x, y, z], dtype=int))

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

            if espressomd.has_features(["PARTICLE_ANISOTROPY"]):
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
        def test_yy_vs_relative(self):
            self.system.part.add(id=0, pos=(0, 0, 0))
            self.system.part.add(id=1, pos=(0, 0, 0))
            self.system.part[1].vs_relative = (0, 5.0, (0.5, -0.5, -0.5, -0.5))
            self.system.part[1].vs_quat = [1, 2, 3, 4]
            np.testing.assert_array_equal(self.system.part[1].vs_quat, [1, 2, 3, 4])
            res = self.system.part[1].vs_relative
            self.assertEqual(res[0], 0, "vs_relative: " + res.__str__())
            self.assertEqual(res[1], 5.0, "vs_relative: " + res.__str__())
            np.testing.assert_allclose(
                res[2], np.array(
                    (0.5, -0.5, -0.5, -0.5)), err_msg="vs_relative: " + res.__str__(), atol=self.tol)

    @ut.skipIf(not espressomd.has_features("DIPOLES"),
               "Features not available, skipping test!")
    def test_contradicting_properties_dip_dipm(self):
        with self.assertRaises(ValueError):
            self.system.part.add(pos=[0, 0, 0], dip=[1, 1, 1], dipm=1.0)

    @ut.skipIf(not espressomd.has_features("DIPOLES", "ROTATION"),
               "Features not available, skipping test!")
    def test_contradicting_properties_dip_quat(self):
        with self.assertRaises(ValueError):
            self.system.part.add(pos=[0, 0, 0], dip=[
                                 1, 1, 1], quat=[1.0, 1.0, 1.0, 1.0])

    @ut.skipIf(not espressomd.has_features("ELECTROSTATICS"), "Test needs ELECTROSTATICS")
    def test_particle_selection(self):
        s = self.system
        s.part.clear()
        positions = ((0.2, 0.3, 0.4), (0.4, 0.2, 0.3), (0.7, 0.7, 0.7))
        charges = 0, 1E-6, -1, 1

        # Place particles
        i = 0
        for pos in positions:
            for q in charges:
                s.part.add(pos=pos, q=q, id=i)
                i += 1

        # Scalar property
        res = s.part.select(q=0)
        self.assertEqual(len(res.id), len(positions))
        for p in res:
            self.assertAlmostEqual(p.q, 0, places=13)

        # Vectorial property
        res = s.part.select(pos=(0.2, 0.3, 0.4))
        self.assertEqual(len(res.id), len(charges))
        for p in res:
            np.testing.assert_allclose(
                (0.2, 0.3, 0.4), np.copy(p.pos), atol=1E-12)

        # Two criteria
        res = s.part.select(pos=(0.2, 0.3, 0.4), q=0)
        self.assertEqual(tuple(res.id), (0,))

        # Emtpy result
        res = s.part.select(q=17)
        self.assertEqual(tuple(res.id), ())
        # User-specified criterion
        res = s.part.select(lambda p: p.pos[0] < 0.5)
        self.assertEqual(tuple(sorted(res.id)), (0, 1, 2, 3, 4, 5, 6, 7))

    def test_image_box(self):
        s = self.system
        s.part.clear()

        pos = 1.5 * s.box_l

        s.part.add(pos=pos)

        np.testing.assert_equal(np.copy(s.part[0].image_box), [1, 1, 1])

    def test_accessing_invalid_id_raises(self):
        self.system.part.clear()
        handle_to_non_existing_particle = self.system.part[42]
        def get_pos_prop(handle):
            return handle.pos
        self.assertRaises(RuntimeError, get_pos_prop, handle_to_non_existing_particle)

    def test_parallel_property_setters(self):
        s= self.system
        s.part.clear()
        s.part.add(pos=s.box_l*np.random.random((100,3)))

        # Copy individual properties of particle 0
        print("If this test hangs, there is an mpi deadlock in a particle property setter." )
        for p in espressomd.particle_data.particle_attributes:
            # Uncomment to identify guilty property
            #print( p)
            
            if not hasattr(s.part[0],p):
                raise Exception("Inconsistency between ParticleHandle and particle_data.particle_attributes")
            try: 
                setattr(s.part[:],p,getattr(s.part[0],p))
            except AttributeError:
                print("Skipping read-only",p)
            # Cause a differtn mpi callback to uncover deadlock immediately
            x=getattr(s.part[:],p)
            

if __name__ == "__main__":
    #print("Features: ", espressomd.features())
    ut.main()
