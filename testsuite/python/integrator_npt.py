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
import numpy as np
import espressomd
import espressomd.integrate


@utx.skipIfMissingFeatures(["NPT", "LENNARD_JONES"])
class IntegratorNPT(ut.TestCase):

    """This tests the NpT integrator interface."""
    system = espressomd.System(box_l=[1.0, 1.0, 1.0])

    def setUp(self):
        self.system.box_l = [5] * 3
        self.system.time_step = 0.01
        self.system.cell_system.skin = 0.25

    def tearDown(self):
        self.system.part.clear()
        self.system.actors.clear()
        self.system.thermostat.turn_off()
        self.system.integrator.set_vv()

    def test_integrator_exceptions(self):
        # invalid parameters should throw exceptions
        with self.assertRaises(RuntimeError):
            self.system.integrator.set_isotropic_npt(ext_pressure=-1, piston=1)
        with self.assertRaises(RuntimeError):
            self.system.integrator.set_isotropic_npt(ext_pressure=1, piston=0)
        with self.assertRaises(RuntimeError):
            self.system.integrator.set_isotropic_npt(ext_pressure=1, piston=-1)
        with self.assertRaises(RuntimeError):
            self.system.integrator.set_isotropic_npt(
                ext_pressure=1, piston=1, direction=[False, False, False])
        with self.assertRaises(Exception):
            self.system.integrator.set_isotropic_npt(
                ext_pressure=1, piston=1, direction=[True, False])

    def test_integrator_recovery(self):
        # the system is still in a valid state after a failure
        system = self.system
        np.random.seed(42)
        npt_params = {'ext_pressure': 0.001, 'piston': 0.001}
        system.box_l = [6] * 3
        system.part.add(pos=np.random.uniform(0, system.box_l[0], (11, 3)))
        system.non_bonded_inter[0, 0].lennard_jones.set_params(
            epsilon=1, sigma=1, cutoff=2**(1 / 6), shift=0.25)
        system.thermostat.set_npt(kT=1.0, gamma0=2, gammav=0.04, seed=42)
        system.integrator.set_isotropic_npt(**npt_params)

        # get the equilibrium box length for the chosen NpT parameters
        system.integrator.run(2000)
        box_l_ref = system.box_l[0]

        # resetting the NpT integrator with incorrect values doesn't leave the
        # system in an undefined state (the old parameters aren't overwritten)
        with self.assertRaises(RuntimeError):
            system.integrator.set_isotropic_npt(ext_pressure=-1, piston=100)
        with self.assertRaises(RuntimeError):
            system.integrator.set_isotropic_npt(ext_pressure=100, piston=-1)
        # the core state is unchanged
        system.integrator.run(500)
        self.assertAlmostEqual(system.box_l[0], box_l_ref, delta=0.15)

        # setting another integrator with incorrect values doesn't leave the
        # system in an undefined state (the old integrator is still active)
        with self.assertRaises(RuntimeError):
            system.integrator.set_steepest_descent(
                f_max=-10, gamma=0, max_displacement=0.1)
        # the interface state is unchanged
        integrator_state = system.integrator.get_state()
        self.assertIsInstance(integrator_state['integrator'],
                              espressomd.integrate.VelocityVerletIsotropicNPT)
        params = integrator_state['integrator'].get_params()
        self.assertEqual(params['ext_pressure'], npt_params['ext_pressure'])
        self.assertEqual(params['piston'], npt_params['piston'])
        # the core state is unchanged
        system.integrator.run(500)
        self.assertAlmostEqual(system.box_l[0], box_l_ref, delta=0.15)

        # setting the NpT integrator with incorrect values doesn't leave the
        # system in an undefined state (the old integrator is still active)
        system.thermostat.turn_off()
        system.integrator.set_vv()
        system.part.clear()
        system.box_l = [5] * 3
        positions_start = np.array([[0, 0, 0], [1., 0, 0]])
        system.part.add(pos=positions_start)
        with self.assertRaises(RuntimeError):
            system.integrator.set_isotropic_npt(ext_pressure=-1, piston=100)
        # the interface state is unchanged
        self.assertIsInstance(system.integrator.get_state()['integrator'],
                              espressomd.integrate.VelocityVerlet)
        # the core state is unchanged
        system.integrator.run(1)
        np.testing.assert_allclose(
            np.copy(system.part.all().pos),
            positions_start + np.array([[-1.2e-3, 0, 0], [1.2e-3, 0, 0]]))

    def run_with_p3m(self, p3m, **npt_kwargs):
        system = self.system
        np.random.seed(42)
        # set up particles
        system.box_l = [6] * 3
        partcl = system.part.add(
            pos=np.random.uniform(
                0, system.box_l[0], (11, 3)))
        if espressomd.has_features("P3M"):
            partcl.q = np.sign(np.arange(-5, 6))
        if espressomd.has_features("DP3M"):
            partcl.dip = tests_common.random_dipoles(11)
        system.non_bonded_inter[0, 0].lennard_jones.set_params(
            epsilon=1, sigma=1, cutoff=2**(1 / 6), shift=0.25)
        system.integrator.set_steepest_descent(
            f_max=10, gamma=0.1, max_displacement=0.01)
        system.integrator.run(100)
        system.integrator.set_vv()
        # combine NpT with a P3M algorithm
        system.actors.add(p3m)
        system.integrator.run(20)
        system.integrator.set_isotropic_npt(ext_pressure=0.001, piston=0.001,
                                            **npt_kwargs)
        system.thermostat.set_npt(kT=1.0, gamma0=2, gammav=0.04, seed=42)
        system.integrator.run(20)

    @utx.skipIfMissingFeatures(["DP3M"])
    def test_dp3m_exception(self):
        # NpT is compatible with DP3M CPU (only cubic box)
        import espressomd.magnetostatics
        dp3m = espressomd.magnetostatics.DipolarP3M(
            prefactor=1.0, accuracy=1e-2, mesh=3 * [36], cao=7, r_cut=1.0,
            alpha=2.995, tune=False)
        with self.assertRaisesRegex(RuntimeError, 'If magnetostatics is being used you must use the cubic box NpT'):
            self.run_with_p3m(
                dp3m, cubic_box=False, direction=(False, True, True))
        self.tearDown()
        try:
            self.run_with_p3m(dp3m)
        except Exception as err:
            self.fail(f'integrator raised ValueError("{err}")')

    @utx.skipIfMissingFeatures(["P3M"])
    def test_p3m_exception(self):
        # NpT is compatible with P3M CPU (only cubic box)
        import espressomd.electrostatics
        p3m = espressomd.electrostatics.P3M(
            prefactor=1.0, accuracy=1e-2, mesh=3 * [8], cao=3, r_cut=0.36,
            alpha=5.35, tune=False)
        with self.assertRaisesRegex(RuntimeError, 'If electrostatics is being used you must use the cubic box NpT'):
            self.run_with_p3m(
                p3m, cubic_box=False, direction=(False, True, True))
        self.tearDown()
        try:
            self.run_with_p3m(p3m)
        except Exception as err:
            self.fail(f'integrator raised ValueError("{err}")')

    @utx.skipIfMissingGPU()
    @utx.skipIfMissingFeatures(["P3M"])
    def test_p3mgpu_exception(self):
        # NpT is not compatible with P3M GPU (no energies)
        import espressomd.electrostatics
        p3m = espressomd.electrostatics.P3MGPU(
            prefactor=1.0, accuracy=1e-2, mesh=3 * [24], cao=2, r_cut=0.24,
            alpha=8.26, tune=False)
        with self.assertRaisesRegex(RuntimeError, 'NpT virial cannot be calculated on P3M GPU'):
            self.run_with_p3m(p3m)


if __name__ == "__main__":
    ut.main()
