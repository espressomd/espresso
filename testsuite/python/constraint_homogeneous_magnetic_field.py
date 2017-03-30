from __future__ import division, print_function

import unittest as ut

import espressomd
import numpy

from espressomd import interactions

@ut.skipIf(not set(["CONSTRAINTS", "DIPOLES"]) < set(espressomd.features()),
        "Features not available, skipping test!")
class HomogeneousMagneticFieldTest(ut.TestCase):

    def test_setter_and_getter(self):
        H_field1 = [0., 1. ,0.]
        H_field2 = [3.533, 5.842, 0.127]

        H_constraint = espressomd.constraints.HomogeneousMagneticField(H=H_field1)
        self.assertEqual(H_constraint.H, H_field1)

        H_constraint.H = H_field2
        for i in range(3):
            self.assertAlmostEqual(H_constraint.H[i], H_field2[i])

    def test_default_value(self):
        H_field_default = [1., 0., 0.] # defined in C++ core
        H_constraint = espressomd.constraints.HomogeneousMagneticField()
        self.assertEqual(H_constraint.H, H_field_default)

    @ut.skipIf("ROTATION" not in espressomd.features(),
            "Feature ROTATION not available, skipping test!")
    def test_add_forces(self):
        S = espressomd.System()
        S.box_l = [3., 3., 3.]
        S.time_step = 0.01
        S.cell_system.skin = 0.4

        H_field = [5., 3., 2.]
        dip_mom = [2., 6., 1.]
        H_constraint = espressomd.constraints.HomogeneousMagneticField(H=H_field)
        S.constraints.add(H_constraint)
        S.part.add(pos=[0,0,0], dip=dip_mom)

        # check that dipole moment is set and no torque is applied
        for i in range(3):
            self.assertEqual(S.part[0].dip[i], dip_mom[i])
            self.assertEqual(S.part[0].torque_lab[i], 0.)

        # check that running the integrator leads to expected torque
        S.integrator.run(0)
        torque_expected = numpy.cross(dip_mom, H_field)
        for i in range(3):
            self.assertAlmostEqual(S.part[0].torque_lab[i], torque_expected[i], places=10)


if __name__ == "__main__":
    ut.main()
