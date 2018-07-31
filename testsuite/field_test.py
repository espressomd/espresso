from __future__ import print_function

import unittest as ut
import numpy as np

import espressomd
from espressomd import constraints

class FieldTest(ut.TestCase):
    """Tests for not space dependend external fields.
    """
    system = espressomd.System(box_l=[1,1,1], time_step=0.01)
    system.cell_system.skin = 0.

    def tearDown(self):
        self.system.constraints.clear()
        self.system.part.clear()

    def test_gravity(self):
        g_const = [1,2,3]
        gravity = constraints.Gravity(g=g_const)

        self.assertSequenceEqual(gravity.g, g_const)
        self.system.constraints.add(gravity)

        if espressomd.has_features("MASS"):
            p = self.system.part.add(pos=[0,0,0], mass=3.1)
        else:
            p = self.system.part.add(pos=[0,0,0])

        self.system.integrator.run(0)

        np.testing.assert_almost_equal(g_const, p.f / p.mass)
        self.assertAlmostEqual(self.system.analysis.energy()['total'], 0.)

    def test_linear_electric_potential(self):
        E = np.array([1.,2.,3.])
        phi0 = 4.

        electric_field = constraints.LinearElectricPotential(E=E, phi0=phi0)
        np.testing.assert_almost_equal(E, electric_field.E)
        self.assertEqual(phi0, electric_field.phi0)

        self.system.constraints.add(electric_field)

        p = self.system.part.add(pos=[0.5,0.5,0.5])

        if espressomd.has_features("ELECTROSTATICS"):
            q_part = -3.1
            p.q = q_part
        else:
            q_part = 0.0

        self.system.integrator.run(0)
        np.testing.assert_almost_equal(q_part * E, p.f)

        self.assertAlmostEqual(self.system.analysis.energy()['total'],
                               q_part * (- np.dot(E, p.pos) + phi0))
        self.assertAlmostEqual(self.system.analysis.energy()['total'],
                               self.system.analysis.energy()['external_fields']

    def test_homogeneous_flow_field(self):
        u = np.array([1.,2.,3.])
        gamma=2.3

        flow_field = constraints.HomogeneousFlowField(u=u, gamma=gamma)
        np.testing.assert_almost_equal(u, flow_field.u)

        self.system.constraints.add(flow_field)

        p = self.system.part.add(pos=[0.5,0.5,0.5], v=[3.,4.,5.])

        self.system.integrator.run(0)
        np.testing.assert_almost_equal(gamma * (u - p.v) , p.f)

        self.assertAlmostEqual(self.system.analysis.energy()['total'],
                               self.system.analysis.energy()['kinetic'])

if __name__ == "__main__":
    print("Features: ", espressomd.features())
    ut.main()

