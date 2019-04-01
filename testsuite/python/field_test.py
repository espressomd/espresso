# Copyright (C) 2010-2018 The ESPResSo project
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
from __future__ import print_function
from itertools import product

import unittest as ut
import numpy as np

import espressomd
from espressomd import constraints


class FieldTest(ut.TestCase):

    """Tests for not space dependend external fields.
    """
    system = espressomd.System(box_l=[10, 10, 10], time_step=0.01)
    system.cell_system.skin = 0.

    def potential(self, x):
        x0 = 5.0 * np.ones_like(x)
        return 0.1 * np.sum(np.power((x - x0), 2))

    def force(self, x):
        x0 = 5.0 * np.ones_like(x)
        return -0.2 * (x - x0)

    def tearDown(self):
        self.system.constraints.clear()
        self.system.part.clear()

    def test_gravity(self):
        g_const = np.array([1, 2, 3])
        gravity = constraints.Gravity(g=g_const)

        np.testing.assert_almost_equal(g_const, np.copy(gravity.g))

        self.system.constraints.add(gravity)

        if espressomd.has_features("MASS"):
            p = self.system.part.add(pos=[0, 0, 0], mass=3.1)
        else:
            p = self.system.part.add(pos=[0, 0, 0])

        self.system.integrator.run(0)

        np.testing.assert_almost_equal(g_const, np.copy(p.f) / p.mass)
        self.assertAlmostEqual(self.system.analysis.energy()['total'], 0.)

    @ut.skipIf(not espressomd.has_features("ELECTROSTATICS"), "Skipping")
    def test_linear_electric_potential(self):
        E = np.array([1., 2., 3.])
        phi0 = 4.

        electric_field = constraints.LinearElectricPotential(E=E, phi0=phi0)
        np.testing.assert_almost_equal(E, electric_field.E)
        self.assertEqual(phi0, electric_field.phi0)

        self.system.constraints.add(electric_field)

        p = self.system.part.add(pos=[0.5, 0.5, 0.5])

        if espressomd.has_features("ELECTROSTATICS"):
            q_part = -3.1
            p.q = q_part
        else:
            q_part = 0.0

        self.system.integrator.run(0)
        np.testing.assert_almost_equal(q_part * E, np.copy(p.f))

        self.assertAlmostEqual(self.system.analysis.energy()['total'],
                               q_part * (- np.dot(E, p.pos) + phi0))
        self.assertAlmostEqual(self.system.analysis.energy()['total'],
                               self.system.analysis.energy()['external_fields'])

    @ut.skipIf(not espressomd.has_features("ELECTROSTATICS"), "Skipping")
    def test_electric_plane_wave(self):
        E0 = np.array([1., -2., 3.])
        k = np.array([-.1, .2, 0.3])
        omega = 5.
        phi = 1.4

        electric_wave = constraints.ElectricPlaneWave(
            E0=E0, k=k, omega=omega, phi=phi)
        np.testing.assert_almost_equal(E0, electric_wave.E0)
        np.testing.assert_almost_equal(k, electric_wave.k)
        np.testing.assert_almost_equal(omega, electric_wave.omega)
        np.testing.assert_almost_equal(phi, electric_wave.phi)

        self.system.constraints.add(electric_wave)

        p = self.system.part.add(pos=[0.4, 0.1, 0.11], q=-14.)
        self.system.time = 1042.

        self.system.integrator.run(0)

        np.testing.assert_almost_equal(np.copy(p.f),
                                       p.q * E0 * np.sin(np.dot(k, p.pos_folded) - omega * self.system.time + phi))

        self.system.integrator.run(10)

        np.testing.assert_almost_equal(np.copy(p.f),
                                       p.q * E0 * np.sin(np.dot(k, p.pos_folded) - omega * self.system.time + phi))

    def test_homogeneous_flow_field(self):
        u = np.array([1., 2., 3.])
        gamma = 2.3

        flow_field = constraints.HomogeneousFlowField(u=u, gamma=gamma)
        np.testing.assert_almost_equal(u, np.copy(flow_field.u))

        self.system.constraints.add(flow_field)

        p = self.system.part.add(pos=[0.5, 0.5, 0.5], v=[3., 4., 5.])

        self.system.integrator.run(0)
        np.testing.assert_almost_equal(gamma * (u - p.v), np.copy(p.f))

        self.assertAlmostEqual(self.system.analysis.energy()['total'],
                               self.system.analysis.energy()['kinetic'])

    def test_potential_field(self):
        h = np.array([.2, .2, .2])
        box = np.array([10., 10., 10.])
        scaling = 2.6

        field_data = constraints.PotentialField.field_from_fn(
            box, h, self.potential)

        F = constraints.PotentialField(field=field_data,
                                       grid_spacing=h,
                                       default_scale=scaling)

        p = self.system.part.add(pos=[0, 0, 0])
        self.system.constraints.add(F)

        for i in product(*map(range, 3 * [10])):
            x = (h * i)
            f_val = F.call_method("_eval_field", x=x)
            np.testing.assert_allclose(f_val, self.potential(x), rtol=1e-3)
            p.pos = x

            self.system.integrator.run(0)
            self.assertAlmostEqual(
                self.system.analysis.energy()['total'], scaling * f_val, places=5)
            np.testing.assert_allclose(
                np.copy(p.f), scaling * self.force(x), rtol=1e-5)

    @ut.skipIf(not espressomd.has_features("ELECTROSTATICS"), "Skipping")
    def test_electric_potential_field(self):
        h = np.array([.2, .2, .2])
        box = np.array([10., 10., 10.])

        field_data = constraints.ElectricPotential.field_from_fn(
            box, h, self.potential)

        F = constraints.ElectricPotential(field=field_data,
                                          grid_spacing=h)

        p = self.system.part.add(pos=[0, 0, 0])
        if espressomd.has_features("ELECTROSTATICS"):
            q_part = -3.1
            p.q = q_part
        else:
            q_part = 0.0

        self.system.constraints.add(F)

        for i in product(*map(range, 3 * [10])):
            x = (h * i)
            f_val = F.call_method("_eval_field", x=x)
            np.testing.assert_allclose(f_val, self.potential(x), rtol=1e-3)
            p.pos = x

            self.system.integrator.run(0)
            self.assertAlmostEqual(
                self.system.analysis.energy()['total'], q_part * f_val, places=5)
            np.testing.assert_allclose(
                np.copy(p.f), q_part * self.force(x), rtol=1e-5)

    def test_force_field(self):
        h = np.array([.8, .8, .8])
        box = np.array([10., 10., 10.])
        scaling = 2.6

        field_data = constraints.ForceField.field_from_fn(box, h, self.force)

        F = constraints.ForceField(field=field_data,
                                   grid_spacing=h,
                                   default_scale=scaling)

        p = self.system.part.add(pos=[0, 0, 0])
        self.system.constraints.add(F)

        for i in product(*map(range, 3 * [10])):
            x = (h * i)
            f_val = np.array(F.call_method("_eval_field", x=x))
            np.testing.assert_allclose(f_val, self.force(x))

            p.pos = x
            self.system.integrator.run(0)
            np.testing.assert_allclose(scaling * f_val, np.copy(p.f))

    def test_flow_field(self):
        h = np.array([.8, .8, .8])
        box = np.array([10., 10., 10.])
        gamma = 2.6

        field_data = constraints.FlowField.field_from_fn(box, h, self.force)

        F = constraints.FlowField(field=field_data,
                                  grid_spacing=h,
                                  gamma=gamma)

        p = self.system.part.add(pos=[0, 0, 0], v=[1, 2, 3])
        self.system.constraints.add(F)

        for i in product(*map(range, 3 * [10])):
            x = (h * i)
            f_val = np.array(F.call_method("_eval_field", x=x))
            np.testing.assert_allclose(f_val, self.force(x))

            p.pos = x
            self.system.integrator.run(0)
            np.testing.assert_allclose(
                -gamma * (p.v - f_val), np.copy(p.f), atol=1e-12)


if __name__ == "__main__":
    print("Features: ", espressomd.features())
    ut.main()
