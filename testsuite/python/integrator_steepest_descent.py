# Copyright (C) 2010-2019 The ESPResSo project
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
import unittest as ut
import unittest_decorators as utx
import numpy as np

import espressomd
import espressomd.constraints


@utx.skipIfMissingFeatures("LENNARD_JONES")
class IntegratorSteepestDescent(ut.TestCase):

    np.random.seed(42)
    system = espressomd.System(box_l=[10.0, 10.0, 10.0])

    test_rotation = espressomd.has_features(("ROTATION", "DIPOLES"))

    box_l = 10.0
    density = 0.6
    vol = box_l**3
    n_part = int(vol * density)

    lj_eps = 1.0
    lj_sig = 1.0
    lj_cut = 2**(1 / 6)

    def setUp(self):
        self.system.box_l = 3 * [self.box_l]
        self.system.cell_system.skin = 0.4
        self.system.time_step = 0.01
        self.system.non_bonded_inter[0, 0].lennard_jones.set_params(
            epsilon=self.lj_eps, sigma=self.lj_sig,
            cutoff=self.lj_cut, shift="auto")
        if self.test_rotation:
            self.system.constraints.add(
                espressomd.constraints.HomogeneousMagneticField(H=[-0.5, 0, 0]))

    def tearDown(self):
        self.system.part.clear()
        self.system.integrator.set_vv()

    def test_relaxation_integrator(self):
        partcls = self.system.part.add(
            pos=np.random.random((self.n_part, 3)) * self.system.box_l)
        if self.test_rotation:
            partcls.dip = np.random.random((self.n_part, 3))
            partcls.dipm = 1
            partcls.rotation = 3 * [True]

        self.assertNotAlmostEqual(
            self.system.analysis.energy()["total"], 0, places=10)

        sd_params = {"f_max": 0.0, "gamma": 0.1, "max_displacement": 0.05}
        self.system.integrator.set_steepest_descent(**sd_params)
        self.system.integrator.run(500)

        self.system.constraints.clear()

        # Check
        self.assertAlmostEqual(
            self.system.analysis.energy()["total"], 0, places=10)
        np.testing.assert_allclose(np.copy(partcls.f), 0.)
        if self.test_rotation:
            np.testing.assert_allclose(np.copy(partcls.dip),
                                       self.n_part * [(-1, 0, 0)], atol=1E-9)

    def test_integration(self):
        max_disp = 0.05
        self.system.part.add(pos=[0, 0, 0], type=0)
        self.system.part.add(pos=[0, 0, self.lj_cut - max_disp / 2], type=0)
        sd_params = {"f_max": 1e-6, "gamma": 0.1, "max_displacement": max_disp}
        self.system.integrator.set_steepest_descent(**sd_params)
        # no displacement for 0 steps
        partcls = self.system.part.all()
        positions = np.copy(partcls.pos)
        steps = self.system.integrator.run(0)
        np.testing.assert_allclose(np.copy(partcls.pos), positions)
        np.testing.assert_allclose(np.copy(partcls.v), 0.)
        np.testing.assert_allclose(partcls.f[:, 0:2], 0.)
        self.assertAlmostEqual(np.sum(partcls.f[:, 2]), 0.)
        self.assertLess(partcls.f[0, 2], -1.)
        self.assertEqual(steps, 0)
        # displacement = max_disp for 1 step
        positions[:, 2] += [-max_disp, max_disp]
        steps = self.system.integrator.run(1)
        np.testing.assert_allclose(np.copy(partcls.pos), positions)
        np.testing.assert_allclose(np.copy(partcls.v), 0.)
        np.testing.assert_allclose(np.copy(partcls.f), 0.)
        self.assertEqual(steps, 1)
        # no displacement after convergence
        steps = self.system.integrator.run(1)
        np.testing.assert_allclose(np.copy(partcls.pos), positions)
        np.testing.assert_allclose(np.copy(partcls.v), 0.)
        np.testing.assert_allclose(np.copy(partcls.f), 0.)
        self.assertEqual(steps, 0)
        # never converges, even when the system is in an energy minimum,
        # because f_max = 0.
        sd_params["f_max"] = 0.0
        self.system.integrator.set_steepest_descent(**sd_params)
        steps = self.system.integrator.run(10)
        np.testing.assert_allclose(np.copy(partcls.pos), positions)
        np.testing.assert_allclose(np.copy(partcls.v), 0.)
        np.testing.assert_allclose(np.copy(partcls.f), 0.)
        self.assertEqual(steps, 10)

    def test_rescaling(self):
        self.system.part.add(pos=[5., 5., 4.9], type=0)
        self.system.part.add(pos=[5., 5., 5.1], type=0)
        partcls = self.system.part.all()

        self.system.non_bonded_inter[0, 0].lennard_jones.set_params(
            epsilon=self.lj_eps, sigma=self.lj_sig,
            cutoff=self.lj_cut, shift="auto")

        self.system.integrator.run(0)
        f_old = np.copy(partcls.f)

        # No-op, because gamma = 0.
        sd_params = {"f_max": 0.0, "gamma": 0.0, "max_displacement": 0.001}

        self.system.integrator.set_steepest_descent(**sd_params)
        self.system.integrator.run(1)

        np.testing.assert_allclose(f_old, np.copy(partcls.f))

    def test_integrator_exceptions(self):
        # invalid parameters should throw exceptions
        with self.assertRaises(RuntimeError):
            self.system.integrator.set_steepest_descent(
                f_max=-1, gamma=1, max_displacement=1)
        with self.assertRaises(RuntimeError):
            self.system.integrator.set_steepest_descent(
                f_max=0, gamma=-1, max_displacement=1)
        with self.assertRaises(RuntimeError):
            self.system.integrator.set_steepest_descent(
                f_max=0, gamma=1, max_displacement=-1)

    def test_integrator_recovery(self):
        # the system is still in a valid state after a failure
        system = self.system
        sd_params = {"f_max": 0, "gamma": 1, "max_displacement": 0.01}
        positions_start = np.array([[0, 0, 0], [1., 0, 0]])
        positions_ref = np.array([[-0.01, 0, 0], [1.01, 0, 0]])
        partcls = system.part.add(pos=positions_start)
        system.integrator.set_steepest_descent(**sd_params)

        # get the positions after one step with the chosen parameters
        system.integrator.run(1)
        positions_ref = np.copy(partcls.pos)

        # resetting the SD integrator with incorrect values doesn't leave the
        # system in an undefined state (the old parameters aren't overwritten)
        with self.assertRaises(RuntimeError):
            system.integrator.set_steepest_descent(
                f_max=-10, gamma=1, max_displacement=0.01)
        with self.assertRaises(RuntimeError):
            system.integrator.set_steepest_descent(
                f_max=0, gamma=-1, max_displacement=0.01)
        with self.assertRaises(RuntimeError):
            system.integrator.set_steepest_descent(
                f_max=0, gamma=1, max_displacement=-1)
        # the core state is unchanged
        partcls.pos = positions_start
        system.integrator.run(1)
        np.testing.assert_allclose(np.copy(partcls.pos), positions_ref)

        # setting another integrator with incorrect values doesn't leave the
        # system in an undefined state (the old integrator is still active)
        if espressomd.has_features("NPT"):
            with self.assertRaises(RuntimeError):
                system.integrator.set_isotropic_npt(ext_pressure=1, piston=-1)
            # the interface state is unchanged
            state = system.integrator.get_state()
            self.assertIsInstance(state['integrator'],
                                  espressomd.integrate.SteepestDescent)
            params = state['integrator'].get_params()
            self.assertEqual(params["f_max"], sd_params["f_max"])
            self.assertEqual(params["gamma"], sd_params["gamma"])
            self.assertEqual(
                params["max_displacement"],
                sd_params["max_displacement"])
            # the core state is unchanged
            partcls.pos = positions_start
            system.integrator.run(1)
            np.testing.assert_allclose(
                np.copy(partcls.pos), positions_ref)

        # setting the SD integrator with incorrect values doesn't leave the
        # system in an undefined state (the old integrator is still active)
        system.integrator.set_vv()
        partcls.pos = positions_start
        with self.assertRaises(RuntimeError):
            system.integrator.set_steepest_descent(
                f_max=0, gamma=1, max_displacement=-1)
        # the interface state is unchanged
        self.assertIsInstance(system.integrator.get_state()['integrator'],
                              espressomd.integrate.VelocityVerlet)
        # the core state is unchanged
        system.integrator.run(1)
        np.testing.assert_allclose(
            np.copy(partcls.pos),
            positions_start + np.array([[-1.2e-3, 0, 0], [1.2e-3, 0, 0]]))


if __name__ == "__main__":
    ut.main()
