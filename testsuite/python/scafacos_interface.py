#
# Copyright (C) 2020-2022 The ESPResSo project
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

import numpy as np
import unittest as ut
import unittest_decorators as utx
import espressomd
import espressomd.electrostatics
import espressomd.magnetostatics
import espressomd.code_info


@utx.skipIfMissingFeatures(["SCAFACOS"])
class ScafacosInterface(ut.TestCase):

    system = espressomd.System(box_l=3 * [5])
    system.time_step = 0.01
    system.cell_system.skin = 0.5
    system.periodicity = 3 * [True]

    def tearDown(self):
        self.system.part.clear()
        self.system.actors.clear()
        self.system.integrator.set_vv()

    def test_decorator(self):
        decorator = utx.skipIfMissingScafacosMethod("unknown")
        with self.assertRaisesRegex(ut.SkipTest, "ScaFaCoS method 'unknown' not available"):
            decorator(lambda: None)()

    def check_available_methods(self, methods):
        # all methods that are accessible from the ScaFaCoS bridge in ESPResSo
        scafacos_methods = {
            "direct", "ewald", "fmm", "memd", "mmm1d", "mmm2d",
            "p2nfft", "p3m", "pepc", "pp3mg", "vmg", "wolf"}

        self.assertGreaterEqual(len(methods), 1)
        for method in methods:
            self.assertIn(method, scafacos_methods)

    @utx.skipIfMissingFeatures(["SCAFACOS_DIPOLES"])
    @utx.skipIfMissingScafacosMethod("p3m")
    @utx.skipIfMissingScafacosMethod("p2nfft")
    def test_magnetostatics_actor_exceptions(self):
        system = self.system
        with self.assertRaisesRegex(ValueError, "Parameter 'prefactor' must be > 0"):
            espressomd.magnetostatics.Scafacos(
                prefactor=-0.8, method_name="p2nfft",
                method_params={"p2nfft_alpha": "0.37"})
        with self.assertRaisesRegex(RuntimeError, "Dipole particles not implemented for solver method 'p3m'"):
            actor = espressomd.magnetostatics.Scafacos(
                prefactor=1., method_name="p3m", method_params={"p3m_cao": 7})
            system.actors.add(actor)
        self.assertEqual(len(system.actors), 0)

    @utx.skipIfMissingFeatures(["SCAFACOS"])
    @utx.skipIfMissingScafacosMethod("p3m")
    @utx.skipIfMissingScafacosMethod("ewald")
    def test_electrostatics_actor_exceptions(self):
        system = self.system
        with self.assertRaisesRegex(ValueError, "Parameter 'prefactor' must be > 0"):
            espressomd.electrostatics.Scafacos(
                prefactor=-0.8, method_name="p3m", method_params={"p3m_cao": 7})
        with self.assertRaisesRegex(ValueError, "Method 'impossible' is unknown or not compiled in ScaFaCoS"):
            espressomd.electrostatics.Scafacos(
                prefactor=1., method_name="impossible", method_params={"p3m_cao": 7})
        with self.assertRaisesRegex(ValueError, "ScaFaCoS methods require at least 1 parameter"):
            espressomd.electrostatics.Scafacos(
                prefactor=1., method_name="p3m", method_params={})

        # choose a method that doesn't support near-field delegation
        scafacos = espressomd.electrostatics.Scafacos(
            prefactor=1.,
            method_name="ewald",
            method_params={"tolerance_field": 0.1})
        self.system.actors.add(scafacos)
        self.assertFalse(scafacos.call_method("get_near_field_delegation"))
        scafacos.call_method("set_near_field_delegation", delegate=False)
        with self.assertRaisesRegex(RuntimeError, "Method 'ewald' cannot delegate short-range calculation"):
            scafacos.call_method("set_near_field_delegation", delegate=True)
        system.actors.clear()

    @utx.skipIfMissingFeatures(["SCAFACOS"])
    @utx.skipIfMissingScafacosMethod("p3m")
    def test_actor_coulomb(self):
        system = self.system

        actor = espressomd.electrostatics.Scafacos(
            prefactor=0.5,
            method_name="p3m",
            method_params={
                "p3m_r_cut": 1.0,
                "p3m_alpha": 2.799269,
                "p3m_grid": 32,
                "p3m_cao": 7})
        system.actors.add(actor)
        params = actor.get_params()
        self.assertEqual(params["prefactor"], 0.5)
        self.assertEqual(params["method_name"], "p3m")
        self.assertEqual(params["method_params"],
                         {'p3m_cao': 7, 'p3m_r_cut': 1.0,
                          'p3m_grid': 32, 'p3m_alpha': 2.799269})

        # check MD cell reset event
        system.box_l = system.box_l
        system.periodicity = system.periodicity

        # force data array update, no-op since there are no particles
        system.integrator.run(0)

        # check available methods
        available_methods = espressomd.code_info.scafacos_methods()
        methods = actor.get_available_methods()
        self.assertEqual(set(methods), set(available_methods))
        self.check_available_methods(methods)

    @utx.skipIfMissingFeatures(["SCAFACOS"])
    @utx.skipIfMissingScafacosMethod("p3m")
    def test_tuning_r_cut_near_field(self):
        system = self.system

        actor = espressomd.electrostatics.Scafacos(
            prefactor=0.5,
            method_name="p3m",
            method_params={
                "p3m_r_cut": -1.0,
                "p3m_alpha": 2.8,
                "p3m_grid": 32,
                "p3m_cao": 7})
        system.actors.add(actor)
        tuned_r_cut = actor.method_params['p3m_r_cut']
        self.assertGreaterEqual(tuned_r_cut, 0.1)
        self.assertLessEqual(tuned_r_cut, min(self.system.box_l) / 2.)

    @utx.skipIfMissingFeatures(["SCAFACOS"])
    @utx.skipIfMissingScafacosMethod("ewald")
    def test_tuning_exceptions(self):
        system = self.system
        # add particles
        N = 100
        np.random.seed(42)
        system.part.add(pos=np.random.uniform(0., system.box_l[0], (N, 3)),
                        q=np.sign((np.arange(N) % 2) * 2. - 1.))

        # minimize system
        system.non_bonded_inter[0, 0].lennard_jones.set_params(
            epsilon=1.0, sigma=1.0, cutoff=2.**(1.0 / 6.0), shift="auto")
        system.integrator.set_steepest_descent(
            f_max=10.,
            gamma=0.001,
            max_displacement=0.01)
        system.integrator.run(100)

        actor = espressomd.electrostatics.Scafacos(
            prefactor=0.5,
            method_name="ewald",
            method_params={"ewald_r_cut": -1.0, "ewald_alpha": 2.8})

        # Ewald cannot delegate near-field calculation
        with self.assertRaisesRegex(RuntimeError, "Method 'ewald' cannot delegate short-range calculation"):
            actor.call_method('set_near_field_delegation', delegate=True)
        actor.call_method('set_near_field_delegation', delegate=False)

        # attempt to tune r_cut with ScaFaCoS
        with self.assertRaisesRegex(RuntimeError, "r_cut is negative"):
            system.actors.add(actor)
        system.actors.clear()

        tuned_r_cut = actor.method_params['ewald_r_cut']
        self.assertAlmostEqual(tuned_r_cut, -1., delta=1e-12)

    @utx.skipIfMissingFeatures(["SCAFACOS"])
    @utx.skipIfMissingScafacosMethod("ewald")
    def test_tuning_r_cut_ewald(self):
        system = self.system

        actor = espressomd.electrostatics.Scafacos(
            prefactor=0.5,
            method_name="ewald",
            method_params={
                "ewald_alpha": 2.,
                "ewald_maxkmax": 200})

        # let ScaFaCoS tune r_cut
        system.actors.add(actor)

        # cutoff is hidden since we don't delegate near-field
        self.assertNotIn("ewald_r_cut", actor.method_params)

    @utx.skipIfMissingFeatures(["SCAFACOS"])
    @utx.skipIfMissingScafacosMethod("ewald")
    def test_tuning_alpha_ewald(self):
        system = self.system

        actor = espressomd.electrostatics.Scafacos(
            prefactor=0.5,
            method_name="ewald",
            method_params={
                "ewald_r_cut": 0.2,
                "ewald_maxkmax": 200})

        # let ScaFaCoS tune alpha
        system.actors.add(actor)

        tuned_r_cut = actor.method_params['ewald_r_cut']
        self.assertAlmostEqual(tuned_r_cut, 0.2, delta=1e-12)

    @utx.skipIfMissingFeatures(["SCAFACOS_DIPOLES"])
    @utx.skipIfMissingScafacosMethod("p2nfft")
    def test_actor_dipoles(self):
        system = self.system

        method_params_ref = {
            "p2nfft_verbose_tuning": 0,
            "pnfft_N": [32, 32, 32],
            "pnfft_n": [32, 32, 32],
            "pnfft_window_name": "bspline",
            "pnfft_m": 4,
            "p2nfft_ignore_tolerance": 1,
            "pnfft_diff_ik": 0,
            "p2nfft_r_cut": 11,
            "p2nfft_alpha": 0.37,
        }
        method_params = {
            "p2nfft_verbose_tuning": "0",  # should come back as an integer
            "pnfft_N": "32,32,32",  # should come back as a list
            "pnfft_n": [32, 32, 32],
            "pnfft_window_name": "bspline",
            "pnfft_m": 4,
            "p2nfft_ignore_tolerance": [1],  # shouldn't come back as a list
            "pnfft_diff_ik": 0,
            "p2nfft_r_cut": 11,  # should come back as an integer
            "p2nfft_alpha": "0.37" + 30 * "0" + "1",  # discard extra digits
        }
        system.actors.add(espressomd.magnetostatics.Scafacos(
            prefactor=1.2,
            method_name="p2nfft",
            method_params=method_params))
        actor = system.actors[0]
        params = actor.get_params()
        self.assertEqual(params["prefactor"], 1.2)
        self.assertEqual(params["method_name"], "p2nfft")
        self.assertEqual(params["method_params"], method_params_ref)
        self.assertEqual({k: type(v) for k, v in actor.method_params.items()},
                         {k: type(v) for k, v in method_params_ref.items()})

        # check MD cell reset event
        system.box_l = system.box_l
        system.periodicity = system.periodicity

        # force data array update, no-op since there are no particles
        system.integrator.run(0)

        # check available methods
        available_methods = espressomd.code_info.scafacos_methods()
        methods = actor.get_available_methods()
        self.assertEqual(set(methods), set(available_methods))
        self.check_available_methods(methods)

    def p3m_data(self):
        system = self.system

        p3m = espressomd.electrostatics.P3M(
            prefactor=0.5,
            accuracy=5e-4,
            mesh=32,
            cao=7,
            r_cut=1.0)
        system.actors.add(p3m)

        dp3m = espressomd.magnetostatics.DipolarP3M(
            prefactor=1.0,
            accuracy=1e-5,
            cao=7,
            mesh=48,
            r_cut=1.88672,
            epsilon="metallic")
        system.actors.add(dp3m)

        system.integrator.run(0, recalc_forces=True)
        ref_E_coulomb = system.analysis.energy()["coulomb"]
        ref_E_dipoles = system.analysis.energy()["dipolar"]
        ref_forces = np.copy(system.part.all().f)
        ref_torques = np.copy(system.part.all().torque_lab)

        system.actors.clear()

        return (ref_E_coulomb, ref_E_dipoles, ref_forces, ref_torques)

    def fcs_data(self):
        system = self.system

        scafacos_coulomb = espressomd.electrostatics.Scafacos(
            prefactor=0.5,
            method_name="p2nfft",
            method_params={
                "p2nfft_verbose_tuning": 0,
                "pnfft_N": "32,32,32",
                "pnfft_n": "32,32,32",
                "tolerance_field": "5e-4",
                "pnfft_window_name": "bspline",
                "pnfft_m": "4",
                "p2nfft_ignore_tolerance": "1",
                "pnfft_diff_ik": "0",
                "p2nfft_r_cut": "1.0",
                "p2nfft_alpha": "2.92"})
        system.actors.add(scafacos_coulomb)

        scafacos_dipoles = espressomd.magnetostatics.Scafacos(
            prefactor=1.0,
            method_name="p2nfft",
            method_params={
                "p2nfft_verbose_tuning": 0,
                "pnfft_N": "32,32,32",
                "pnfft_n": "32,32,32",
                "pnfft_window_name": "bspline",
                "pnfft_m": "4",
                "p2nfft_ignore_tolerance": "1",
                "pnfft_diff_ik": "0",
                "p2nfft_r_cut": "11",
                "p2nfft_alpha": "0.37"})
        system.actors.add(scafacos_dipoles)

        system.integrator.run(0, recalc_forces=True)
        ref_E_coulomb = system.analysis.energy()["coulomb"]
        ref_E_dipoles = system.analysis.energy()["dipolar"]
        ref_forces = np.copy(system.part.all().f)
        ref_torques = np.copy(system.part.all().torque_lab)

        # check MD cell reset has no impact
        system.box_l = system.box_l
        system.periodicity = system.periodicity
        system.cell_system.node_grid = system.cell_system.node_grid
        system.integrator.run(0, recalc_forces=True)
        new_E_coulomb = system.analysis.energy()["coulomb"]
        new_E_dipoles = system.analysis.energy()["dipolar"]
        new_forces = np.copy(system.part.all().f)
        new_torques = np.copy(system.part.all().torque_lab)
        self.assertAlmostEqual(new_E_coulomb, ref_E_coulomb, delta=0)
        self.assertAlmostEqual(new_E_dipoles, ref_E_dipoles, delta=0)
        np.testing.assert_allclose(new_forces, ref_forces, atol=0, rtol=0.)
        np.testing.assert_allclose(new_torques, ref_torques, atol=0, rtol=0.)

        system.actors.clear()

        return (ref_E_coulomb, ref_E_dipoles, ref_forces, ref_torques)

    @utx.skipIfMissingFeatures(["LENNARD_JONES", "P3M", "SCAFACOS_DIPOLES"])
    @utx.skipIfMissingScafacosMethod("p2nfft")
    def test_electrostatics_plus_magnetostatics(self):
        # check that two instances of ScaFaCoS can be used
        system = self.system

        # add particles
        N = 100
        np.random.seed(42)
        system.part.add(pos=np.random.uniform(0., system.box_l[0], (N, 3)),
                        dip=np.random.uniform(0., 1., (N, 3)),
                        q=np.sign((np.arange(N) % 2) * 2. - 1.),
                        rotation=N * [(True, True, True)])

        # minimize system
        system.non_bonded_inter[0, 0].lennard_jones.set_params(
            epsilon=1.0, sigma=1.0, cutoff=2.**(1.0 / 6.0), shift="auto")
        system.integrator.set_steepest_descent(
            f_max=10.,
            gamma=0.001,
            max_displacement=0.01)
        system.integrator.run(100)
        system.integrator.set_vv()

        # compute forces and energies
        p3m_E_coulomb, p3m_E_dipoles, p3m_forces, p3m_torques = self.p3m_data()
        fcs_E_coulomb, fcs_E_dipoles, fcs_forces, fcs_torques = self.fcs_data()

        self.assertAlmostEqual(fcs_E_coulomb, p3m_E_coulomb, delta=1e-4)
        self.assertAlmostEqual(fcs_E_dipoles, p3m_E_dipoles, delta=1e-4)
        np.testing.assert_allclose(fcs_forces, p3m_forces, rtol=1e-2)
        np.testing.assert_allclose(fcs_torques, p3m_torques, rtol=1e-3)


if __name__ == "__main__":
    ut.main()
