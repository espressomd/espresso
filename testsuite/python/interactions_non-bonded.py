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
import espressomd
import numpy as np
import unittest as ut
import unittest_decorators as utx
import tests_common


class InteractionsNonBondedTest(ut.TestCase):
    system = espressomd.System(box_l=3 * [10.])
    system.cell_system.skin = 0.
    system.time_step = .1

    start_pos = np.random.rand(3) * system.box_l
    axis = np.random.rand(3)
    axis /= np.linalg.norm(axis)
    step = axis * 0.01
    step_width = np.linalg.norm(step)

    def setUp(self):
        self.system.part.add(pos=self.start_pos, type=0)
        self.system.part.add(pos=self.start_pos, type=0)

    def tearDown(self):
        self.system.non_bonded_inter.reset()
        self.system.part.clear()

    # Required, since assertAlmostEqual does NOT check significant places
    def assertFractionAlmostEqual(self, a, b, **args):
        if abs(b) < 1E-8:
            self.assertAlmostEqual(a, b, **args)
        else:
            self.assertAlmostEqual(a / b, 1., **args)

    def assertItemsFractionAlmostEqual(self, a, b):
        for i, ai in enumerate(a):
            self.assertFractionAlmostEqual(ai, b[i])

    #
    # Tests
    #

    # Test Generic Lennard-Jones Potential
    @utx.skipIfMissingFeatures("LENNARD_JONES_GENERIC")
    def test_lj_generic(self):

        self.run_test("generic_lennard_jones",
                      {"epsilon": 2.12,
                       "sigma": 1.37,
                       "cutoff": 2.122,
                       "offset": 0.185,
                       "b1": 4.22,
                       "b2": 3.63,
                       "e1": 10.32,
                       "e2": 5.81,
                       "shift": -0.13},
                      force_kernel=tests_common.lj_generic_force,
                      energy_kernel=tests_common.lj_generic_potential,
                      n_steps=231,
                      force_kernel_needs_espressomd=True)

    # Test WCA Potential
    @utx.skipIfMissingFeatures("WCA")
    def test_wca(self):

        wca_eps = 2.12
        wca_sig = 1.37
        wca_cutoff = wca_sig * 2.**(1. / 6.)
        wca_shift = -((wca_sig / wca_cutoff)**12 - (wca_sig / wca_cutoff)**6)

        self.run_test("wca",
                      {"epsilon": wca_eps,
                       "sigma": wca_sig},
                      force_kernel=lambda espressomd, r, epsilon, sigma: tests_common.lj_generic_force(
                          espressomd, r, epsilon=epsilon, sigma=sigma, cutoff=wca_cutoff),
                      energy_kernel=lambda r, epsilon, sigma: tests_common.lj_generic_potential(
                          r, epsilon=epsilon, sigma=sigma, cutoff=wca_cutoff, shift=4. * wca_shift),
                      n_steps=231,
                      force_kernel_needs_espressomd=True)

    # Test Generic Lennard-Jones Softcore Potential
    @utx.skipIfMissingFeatures("LJGEN_SOFTCORE")
    def test_lj_generic_softcore(self):

        self.run_test("generic_lennard_jones",
                      {"epsilon": 2.12,
                       "sigma": 1.37,
                       "cutoff": 2.125,
                       "offset": 0.182,
                       "b1": 6.22,
                       "b2": 3.63,
                       "e1": 13.32,
                       "e2": 3.74,
                       "shift": 0.13,
                       "delta": 0.1,
                       "lam": 0.34},
                      force_kernel=tests_common.lj_generic_force,
                      energy_kernel=tests_common.lj_generic_potential,
                      n_steps=231,
                      force_kernel_needs_espressomd=True)

    # Test Lennard-Jones Potential
    @utx.skipIfMissingFeatures("LENNARD_JONES")
    def test_lj(self):

        self.run_test("lennard_jones",
                      {"epsilon": 1.92,
                       "sigma": 1.03,
                       "cutoff": 1.123,
                       "shift": 0.92},
                      force_kernel=tests_common.lj_force,
                      energy_kernel=tests_common.lj_potential,
                      n_steps=113,
                      force_kernel_needs_espressomd=True)

    # Test Lennard-Jones Cosine Potential
    @utx.skipIfMissingFeatures("LJCOS")
    def test_lj_cos(self):

        self.run_test("lennard_jones_cos",
                      {"epsilon": 3.32,
                       "sigma": 0.73,
                       "cutoff": 1.523,
                       "offset": 0.223},
                      force_kernel=tests_common.lj_cos_force,
                      energy_kernel=tests_common.lj_cos_potential,
                      n_steps=175,
                      force_kernel_needs_espressomd=True)

    # Test Lennard-Jones Cosine^2 Potential
    @utx.skipIfMissingFeatures("LJCOS2")
    def test_lj_cos2(self):

        self.run_test("lennard_jones_cos2",
                      {"epsilon": 0.31,
                       "sigma": 0.73,
                       "width": 1.523,
                       "offset": 0.321},
                      force_kernel=tests_common.lj_cos2_force,
                      energy_kernel=tests_common.lj_cos2_potential,
                      n_steps=267,
                      force_kernel_needs_espressomd=True)

    # Test Smooth-step Potential
    @utx.skipIfMissingFeatures("SMOOTH_STEP")
    def test_smooth_step(self):

        self.run_test("smooth_step",
                      {"eps": 4.92,
                       "sig": 3.03,
                       "cutoff": 1.253,
                       "d": 2.52,
                       "n": 11,
                       "k0": 2.13},
                      force_kernel=tests_common.smooth_step_force,
                      energy_kernel=tests_common.smooth_step_potential,
                      n_steps=126)

    # Test BMHTF Potential
    @utx.skipIfMissingFeatures("BMHTF_NACL")
    def test_bmhtf(self):

        self.run_test("bmhtf",
                      {"a": 3.92,
                       "b": 2.43,
                       "c": 1.23,
                       "d": 3.33,
                       "sig": 0.123,
                       "cutoff": 1.253},
                      force_kernel=tests_common.bmhtf_force,
                      energy_kernel=tests_common.bmhtf_potential,
                      n_steps=126)

    # Test Morse Potential
    @utx.skipIfMissingFeatures("MORSE")
    def test_morse(self):

        self.run_test("morse",
                      {"eps": 1.92,
                       "alpha": 3.03,
                       "rmin": 0.123,
                       "cutoff": 1.253},
                      force_kernel=tests_common.morse_force,
                      energy_kernel=tests_common.morse_potential,
                      n_steps=126)

    # Test Buckingham Potential
    @utx.skipIfMissingFeatures("BUCKINGHAM")
    def test_buckingham(self):

        self.run_test("buckingham",
                      {"a": 3.71,
                       "b": 2.92,
                       "c": 5.32,
                       "d": 4.11,
                       "discont": 1.03,
                       "cutoff": 2.253,
                       "shift": 0.133},
                      force_kernel=tests_common.buckingham_force,
                      energy_kernel=tests_common.buckingham_potential,
                      n_steps=226,
                      force_kernel_remove_shift=False)

    # Test Soft-sphere Potential
    @utx.skipIfMissingFeatures("SOFT_SPHERE")
    def test_soft_sphere(self):

        self.run_test("soft_sphere",
                      {"a": 1.92,
                       "n": 3.03,
                       "cutoff": 1.123,
                       "offset": 0.123},
                      force_kernel=tests_common.soft_sphere_force,
                      energy_kernel=tests_common.soft_sphere_potential,
                      n_steps=113,
                      n_initial_steps=12)

    # Test Hertzian Potential
    @utx.skipIfMissingFeatures("HERTZIAN")
    def test_hertzian(self):

        self.run_test("hertzian",
                      {"eps": 6.92,
                       "sig": 2.432},
                      force_kernel=tests_common.hertzian_force,
                      energy_kernel=tests_common.hertzian_potential,
                      n_steps=244)

    # Test Gaussian Potential
    @utx.skipIfMissingFeatures("GAUSSIAN")
    def test_gaussian(self):

        self.run_test("gaussian",
                      {"eps": 6.92,
                       "sig": 4.03,
                       "cutoff": 1.243},
                      force_kernel=tests_common.gaussian_force,
                      energy_kernel=tests_common.gaussian_potential,
                      n_steps=125)

    # Test the Gay-Berne potential and the resulting force and torque
    @utx.skipIfMissingFeatures("GAY_BERNE")
    def test_gb(self):

        # helper function definitions
        def gradient(func, x0, dx=1.0e-7):
            """
            Approximate the gradient of a function at a point x0
            using the two-point central difference formula with spacing 2dx.

            Parameters
            ----------
            func: :obj:`function`
                function for which the gradient is calculated
            x0: (3,) array_like of :obj:`float`
                Point in N-dimensional space where the derivatives are calculated
            dx: :obj:`float`, optional
                Spacing
            Returns
            -------
            (3,) array_like of obj:`float`
                the approximated gradient of func at x0

            """
            def partial_x(x):
                return (func(x0 + x) - func(x0 - x)) / (
                    2.0 * np.linalg.norm(x))
            delta = np.array([dx, 0.0, 0.0])
            return np.array([partial_x(np.roll(delta, i)) for i in range(3)])

        def setup_system(gb_params):
            k_1, k_2, mu, nu, sigma_0, epsilon_0, cut = gb_params

            self.system.part.clear()
            self.system.part.add(
                pos=(1, 2, 3), rotation=(1, 1, 1), type=0)
            self.system.part.add(
                pos=(2.2, 2.1, 2.9), rotation=(1, 1, 1), type=0)

            self.system.non_bonded_inter[0, 0].gay_berne.set_params(
                sig=sigma_0, cut=cut, eps=epsilon_0, k1=k_1, k2=k_2, mu=mu,
                nu=nu)

        def advance_and_rotate_part(particle):
            particle.pos = particle.pos + self.step
            particle.rotate(axis=(1, 2, 3), angle=0.3)
            particle.rotate(axis=(1, -2, -4), angle=1.2)

        def get_simulation_energy():
            return self.system.analysis.energy()["non_bonded"]

        def get_reference_energy(gb_params, r, director1, director2):
            k_1, k_2, mu, nu, sigma_0, epsilon_0, cut = gb_params
            r_cut = r * cut / np.linalg.norm(r)

            E_ref = tests_common.gay_berne_potential(
                r, director1, director2, epsilon_0, sigma_0, mu, nu, k_1, k_2)

            E_ref -= tests_common.gay_berne_potential(
                r_cut, director1, director2, epsilon_0, sigma_0, mu, nu,
                k_1, k_2)
            return E_ref

        def get_reference_force(gb_params, r, dir1, dir2):
            return -gradient(
                lambda x: get_reference_energy(gb_params, x, dir1, dir2),
                x0=r, dx=1.0e-7)

        def get_reference_torque(gb_params, r, dir1, dir2):
            force_in_dir1 = gradient(
                lambda x: get_reference_energy(gb_params, r, x, dir2),
                x0=dir1, dx=1.0e-7)

            return np.cross(-dir1, force_in_dir1)

        # actual tests of the gb potential

        k_1 = 1.2
        k_2 = 2.4
        mu = 2.
        nu = 5.
        sigma_0 = 1.2
        epsilon_0 = 0.8
        cut = 3.3

        gb_params = (k_1, k_2, mu, nu, sigma_0, epsilon_0, cut)

        setup_system(gb_params)

        p1, p2 = self.system.part.all()

        delta = 1.0e-6

        for _ in range(100):

            advance_and_rotate_part(p2)
            self.system.integrator.run(recalc_forces=True, steps=0)

            r = self.system.distance_vec(p1, p2)
            director1 = p1.director
            director2 = p2.director

            # Calc energies
            E_sim = get_simulation_energy()
            E_ref = get_reference_energy(gb_params, r, director1, director2)
            # Test energies
            self.assertAlmostEqual(E_sim, E_ref, delta=delta)

            # Calc forces
            f1_sim = np.copy(p1.f)
            f2_sim = np.copy(p2.f)
            f2_ref = get_reference_force(gb_params, r, director1, director2)
            # Test forces
            # force equals minus the counter-force
            np.testing.assert_array_equal(f1_sim, -f2_sim)
            # compare force to reference force
            for i in range(3):
                self.assertAlmostEqual(f2_sim[i], f2_ref[i], delta=delta)

            # Calc torques
            torque1_sim = p1.torque_lab
            torque2_sim = p2.torque_lab
            torque1_ref = get_reference_torque(
                gb_params, r, director1, director2)
            torque2_ref = get_reference_torque(
                gb_params, r, director2, director1)
            # Test torques
            for i in range(3):
                self.assertAlmostEqual(
                    torque1_sim[i],
                    torque1_ref[i],
                    delta=delta)
                self.assertAlmostEqual(
                    torque2_sim[i],
                    torque2_ref[i],
                    delta=delta)

        # Test zero energy
        self.system.non_bonded_inter[0, 0].gay_berne.set_params(
            sig=sigma_0, cut=0, eps=0, k1=k_1, k2=k_2, mu=mu, nu=nu)
        self.system.integrator.run(0)
        self.assertEqual(self.system.analysis.energy()["non_bonded"], 0.0)

    def run_test(self, name, parameters, force_kernel,
                 energy_kernel, n_steps, n_initial_steps=0,
                 force_kernel_needs_espressomd=False,
                 force_kernel_remove_shift=True):

        getattr(self.system.non_bonded_inter[0, 0], name).set_params(
            **parameters)
        p0, p1 = self.system.part.all()
        p1.pos = p0.pos + self.step * n_initial_steps

        force_parameters = parameters.copy()
        if "shift" in force_parameters and force_kernel_remove_shift:
            del force_parameters["shift"]
        if force_kernel_needs_espressomd:
            force_parameters["espressomd"] = espressomd

        for _ in range(n_steps):
            p1.pos = p1.pos + self.step
            d = np.linalg.norm(p1.pos - p0.pos)
            self.system.integrator.run(recalc_forces=True, steps=0)

            # Calculate energies
            E_sim = self.system.analysis.energy()["non_bonded"]
            E_ref = energy_kernel(r=d, **parameters)

            # Calculate forces
            f0_sim = np.copy(p0.f)
            f1_sim = np.copy(p1.f)
            f1_ref = self.axis * force_kernel(r=d, **force_parameters)

            # Check that energies match ...
            self.assertFractionAlmostEqual(E_sim, E_ref)
            # force equals minus the counter-force ...
            np.testing.assert_array_equal(f0_sim, -f1_sim)
            # and has correct value.
            self.assertItemsFractionAlmostEqual(f1_sim, f1_ref)

        # forces and energies are zero beyond the interaction cutoff
        p1.pos = p0.pos + self.system.box_l / 2
        self.system.integrator.run(recalc_forces=True, steps=0)
        E_sim = self.system.analysis.energy()["non_bonded"]
        np.testing.assert_array_equal(E_sim, 0.)
        np.testing.assert_array_equal(np.copy(p0.f), 0.)
        np.testing.assert_array_equal(np.copy(p1.f), 0.)


if __name__ == '__main__':
    ut.main()
