#
# Copyright (C) 2020 The ESPResSo project
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
import espressomd.scafacos


@utx.skipIfMissingFeatures(["SCAFACOS"])
class ScafacosInterface(ut.TestCase):

    system = espressomd.System(box_l=3 * [5])
    system.time_step = 0.01
    system.cell_system.skin = 0.5
    system.periodicity = 3 * [True]

    def tearDown(self):
        self.system.part.clear()
        self.system.actors.clear()

    def test_available_methods(self):
        # all methods that are accessible from the ScaFaCoS bridge in ESPResSo
        scafacos_methods = {
            "direct", "ewald", "fmm", "memd", "mmm1d", "mmm2d",
            "p2nfft", "p3m", "pepc", "pp3mg", "vmg", "wolf"}

        # all methods that were compiled in when building ScaFaCoS
        available_methods = espressomd.scafacos.available_methods()
        self.assertGreaterEqual(len(available_methods), 1)
        for method in available_methods:
            self.assertIn(method, scafacos_methods)

    @ut.skipIf(not espressomd.has_features('SCAFACOS') or
               'p3m' not in espressomd.scafacos.available_methods(),
               'Skipping test: missing ScaFaCoS p3m method')
    def test_actor_exceptions(self):
        system = self.system

        if espressomd.has_features('SCAFACOS_DIPOLES'):
            with self.assertRaisesRegex(ValueError, "Dipolar prefactor has to be >= 0"):
                system.actors.add(espressomd.magnetostatics.Scafacos(
                    prefactor=-1, method_name="p3m", method_params={"p3m_cao": 7}))
            system.actors.clear()
        with self.assertRaisesRegex(ValueError, "Coulomb prefactor has to be >= 0"):
            system.actors.add(espressomd.electrostatics.Scafacos(
                prefactor=-1, method_name="p3m", method_params={"p3m_cao": 7}))
        system.actors.clear()
        with self.assertRaisesRegex(ValueError, "method 'impossible' is unknown or not compiled in ScaFaCoS"):
            system.actors.add(espressomd.electrostatics.Scafacos(
                prefactor=1, method_name="impossible", method_params={"p3m_cao": 7}))
        system.actors.clear()
        with self.assertRaisesRegex(ValueError, "ScaFaCoS methods require at least 1 parameter"):
            system.actors.add(espressomd.electrostatics.Scafacos(
                prefactor=1, method_name="p3m", method_params={}))
        system.actors.clear()

    @ut.skipIf(not espressomd.has_features('SCAFACOS') or
               'p3m' not in espressomd.scafacos.available_methods(),
               'Skipping test: missing ScaFaCoS p3m method')
    def test_actor_coulomb(self):
        system = self.system

        system.actors.add(espressomd.electrostatics.Scafacos(
            prefactor=0.5,
            method_name="p3m",
            method_params={
                "p3m_r_cut": 1.0,
                "p3m_alpha": 2.799269,
                "p3m_grid": 32,
                "p3m_cao": 7}))
        actor = system.actors[0]
        params = actor.get_params()
        self.assertEqual(params["prefactor"], 0.5)
        self.assertEqual(params["method_name"], "p3m")
        self.assertEqual(params["method_params"],
                         {'p3m_cao': '7', 'p3m_r_cut': '1.0',
                          'p3m_grid': '32', 'p3m_alpha': '2.799269'})

        # check MD cell reset event
        system.box_l = system.box_l
        system.periodicity = system.periodicity

    @ut.skipIf(not espressomd.has_features('SCAFACOS_DIPOLES') or
               'p2nfft' not in espressomd.scafacos.available_methods(),
               'Skipping test: missing ScaFaCoS p2nfft method')
    def test_actor_dipoles(self):
        system = self.system

        method_params = {
            "p2nfft_verbose_tuning": "0",
            "pnfft_N": "32,32,32",
            "pnfft_n": "32,32,32",
            "pnfft_window_name": "bspline",
            "pnfft_m": "4",
            "p2nfft_ignore_tolerance": "1",
            "pnfft_diff_ik": "0",
            "p2nfft_r_cut": "11",
            "p2nfft_alpha": "0.37"}
        system.actors.add(espressomd.magnetostatics.Scafacos(
            prefactor=1.2,
            method_name="p2nfft",
            method_params=method_params))
        actor = system.actors[0]
        params = actor.get_params()
        self.assertEqual(params["prefactor"], 1.2)
        self.assertEqual(params["method_name"], "p2nfft")
        self.assertEqual(params["method_params"], method_params)

        # check MD cell reset event
        system.box_l = system.box_l
        system.periodicity = system.periodicity

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

        system.actors.clear()

        return (ref_E_coulomb, ref_E_dipoles, ref_forces, ref_torques)

    @utx.skipIfMissingFeatures(["LENNARD_JONES", "P3M"])
    @ut.skipIf(not espressomd.has_features('SCAFACOS_DIPOLES') or
               'p2nfft' not in espressomd.scafacos.available_methods(),
               'Skipping test: missing SCAFACOS_DIPOLES or p2nfft method')
    def test_electrostatics_plus_magnetostatics(self):
        # check that two instances of ScaFaCoS can be used
        system = self.system

        # add particles
        N = 100
        np.random.seed(42)
        system.part.add(pos=np.random.uniform(0, system.box_l[0], (N, 3)),
                        dip=np.random.uniform(0, 1, (N, 3)),
                        q=np.sign((np.arange(N) % 2) * 2 - 1),
                        rotation=N * [(1, 1, 1)])

        # minimize system
        system.non_bonded_inter[0, 0].lennard_jones.set_params(
            epsilon=1.0, sigma=1.0, cutoff=2**(1.0 / 6.0), shift="auto")
        system.integrator.set_steepest_descent(
            f_max=10,
            gamma=0.001,
            max_displacement=0.01)
        system.integrator.run(100)

        # compute forces and energies
        p3m_E_coulomb, p3m_E_dipoles, p3m_forces, p3m_torques = self.p3m_data()
        fcs_E_coulomb, fcs_E_dipoles, fcs_forces, fcs_torques = self.fcs_data()

        self.assertAlmostEqual(fcs_E_coulomb, p3m_E_coulomb, delta=1e-4)
        self.assertAlmostEqual(fcs_E_dipoles, p3m_E_dipoles, delta=1e-4)
        np.testing.assert_allclose(fcs_forces, p3m_forces, rtol=1e-2)
        np.testing.assert_allclose(fcs_torques, p3m_torques, rtol=1e-3)


if __name__ == "__main__":
    ut.main()
