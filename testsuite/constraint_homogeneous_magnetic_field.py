from __future__ import division, print_function

import unittest as ut

import espressomd
import numpy
import sys

class HomogeneousMagneticFieldTest(ut.TestCase):

    S = espressomd.System(box_l=[1.0, 1.0, 1.0])
    S.seed  = S.cell_system.get_state()['n_nodes'] * [1234]
    numpy.random.seed(S.seed)


    def setUp(self):
        self.S.box_l = [3.0, 3.0, 3.0]
        self.S.time_step = 0.01
        self.S.cell_system.skin = 0.4

    def tearDown(self):
        self.S.constraints.clear()

    def test_setter_and_getter(self):
        H_field1 = [0.0, 1.0, 0.0]
        H_field2 = [3.533, 5.842, 0.127]

        H_constraint = espressomd.constraints.HomogeneousMagneticField(
            H=H_field1)
        self.assertEqual(H_constraint.H, H_field1)

        H_constraint.H = H_field2
        for i in range(3):
            self.assertAlmostEqual(H_constraint.H[i], H_field2[i])

    def test_default_value(self):
        H_field_default = [1.0, 0.0, 0.0]  # defined in C++ core
        H_constraint = espressomd.constraints.HomogeneousMagneticField()
        self.assertEqual(H_constraint.H, H_field_default)

    @ut.skipIf(not espressomd.has_features(["DIPOLES"]),
               "Features DIPOLES not available, skipping test!")
    def test_add_energy_and_forces(self):
        H_field = [5.0, 3.0, 2.0]
        dip_mom0 = [2.0, 6.0, 1.]
        dip_mom1 = [-1.0, 0.5, -0.2]

        # check that the dipolar energy is zero initially, ...
        self.assertEqual(self.S.analysis.energy()["dipolar"], 0.0)

        H_constraint = espressomd.constraints.HomogeneousMagneticField(
            H=H_field)
        self.S.constraints.add(H_constraint)

        # ... and also after adding the constraint
        self.assertEqual(self.S.analysis.energy()["dipolar"], 0.0)

        # check dipolar energy when adding dipole moments
        self.S.part.add(id=0, pos=[0,0,0], dip=dip_mom0,rotation=(1,1,1))
        self.assertEqual(self.S.analysis.energy()["dipolar"], -1.0*numpy.dot(H_field, dip_mom0))
        self.S.part.add(id=1, pos=[1,1,1], dip=dip_mom1,rotation=(1,1,1))
        self.assertEqual(self.S.analysis.energy()["dipolar"],
                         (-1.0 * numpy.dot(H_field, dip_mom0) - 1.0 * numpy.dot(H_field, dip_mom1)))

        if espressomd.has_features(["ROTATION"]):
            # check that running the integrator leads to expected torques
            self.S.integrator.run(0)
            torque_expected0 = numpy.cross(dip_mom0, H_field)
            torque_expected1 = numpy.cross(dip_mom1, H_field)
            for i in range(3):
                self.assertAlmostEqual(
                    self.S.part[0].torque_lab[i],
                    torque_expected0[i],
                    places=10)
                self.assertAlmostEqual(
                    self.S.part[1].torque_lab[i],
                    torque_expected1[i],
                    places=10)


if __name__ == "__main__":
    print("Features: ", espressomd.features())
    ut.main()
