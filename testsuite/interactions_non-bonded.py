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
from __future__ import print_function
import espressomd
import numpy
import unittest as ut
import tests_common


class InteractionsNonBondedTest(ut.TestCase):
    system = espressomd.System(box_l=[1.0, 1.0, 1.0])
    box_l = 10.

    start_pos = numpy.random.rand(3) * box_l
    axis = numpy.random.rand(3)
    axis /= numpy.linalg.norm(axis)
    step = axis * 0.01
    step_width = numpy.linalg.norm(step)

    def setUp(self):

        self.system.box_l = [self.box_l] * 3
        self.system.cell_system.skin = 0.4
        self.system.time_step = .1

        self.system.part.add(id=0, pos=self.start_pos, type=0)
        self.system.part.add(id=1, pos=self.start_pos, type=0)

    def tearDown(self):

        self.system.part.clear()

    # Required, since assertAlmostEqual does NOT check significant places
    def assertFractionAlmostEqual(self, a, b):
        if abs(b) < 1E-8:
            self.assertAlmostEqual(a, b)
        else:
            self.assertAlmostEqual(a / b, 1.)

    def assertItemsFractionAlmostEqual(self, a, b):
        for i, ai in enumerate(a):
            self.assertFractionAlmostEqual(ai, b[i])


    #
    # Tests
    #

    # Test Generic Lennard-Jones Potential
    @ut.skipIf(not espressomd.has_features("LENNARD_JONES_GENERIC"),
               "Features not available, skipping test!")
    def test_lj_generic(self):

        lj_eps = 2.12
        lj_sig = 1.37
        lj_cut = 2.122
        lj_off = 0.185
        lj_b1 = 4.22
        lj_b2 = 3.63
        lj_e1 = 10.32
        lj_e2 = 5.81
        lj_shift = -0.13

        self.system.non_bonded_inter[0, 0].generic_lennard_jones.set_params(
            epsilon=lj_eps, sigma=lj_sig, cutoff=lj_cut, offset=lj_off, b1=lj_b1, b2=lj_b2, e1=lj_e1, e2=lj_e2, shift=lj_shift)

        for i in range(231):
            self.system.part[1].pos = self.system.part[1].pos + self.step
            self.system.integrator.run(recalc_forces=True, steps=0)

            # Calculate energies
            E_sim = self.system.analysis.energy()["non_bonded"]
            E_ref = tests_common.lj_generic_potential(
                r=(i + 1) * self.step_width, eps=lj_eps, sig=lj_sig,
                cutoff=lj_cut, offset=lj_off, b1=lj_b1, b2=lj_b2, e1=lj_e1, e2=lj_e2, shift=lj_shift)

            # Calculate forces
            f0_sim = self.system.part[0].f
            f1_sim = self.system.part[1].f
            f1_ref = self.axis * tests_common.lj_generic_force(espressomd, 
                r=(i + 1) * self.step_width, eps=lj_eps,
                sig=lj_sig, cutoff=lj_cut, offset=lj_off, b1=lj_b1, b2=lj_b2, e1=lj_e1, e2=lj_e2)

            # Check that energies match, ...
            self.assertFractionAlmostEqual(E_sim, E_ref)
            # force equals minus the counter-force  ...
            self.assertTrue((f0_sim == -f1_sim).all())
            # and has correct value.
            self.assertItemsFractionAlmostEqual(f1_sim, f1_ref)

        self.system.non_bonded_inter[0, 0].generic_lennard_jones.set_params(
            epsilon=0.)

    # Test Generic Lennard-Jones Softcore Potential
    @ut.skipIf(not espressomd.has_features("LJGEN_SOFTCORE"),
               "Features not available, skipping test!")
    def test_lj_generic_softcore(self):

        lj_eps = 2.12
        lj_sig = 1.37
        lj_cut = 2.125
        lj_off = 0.182
        lj_b1 = 6.22
        lj_b2 = 3.63
        lj_e1 = 13.32
        lj_e2 = 3.74
        lj_shift = 0.13
        lj_delta = 0.1
        lj_lam = 0.34

        self.system.non_bonded_inter[0, 0].generic_lennard_jones.set_params(
            epsilon=lj_eps, sigma=lj_sig, cutoff=lj_cut, offset=lj_off, b1=lj_b1, b2=lj_b2, e1=lj_e1, e2=lj_e2, shift=lj_shift, delta=lj_delta, lam=lj_lam)

        for i in range(231):
            self.system.part[1].pos = self.system.part[1].pos + self.step
            self.system.integrator.run(recalc_forces=True, steps=0)

            # Calculate energies
            E_sim = self.system.analysis.energy()["non_bonded"]
            E_ref = tests_common.lj_generic_potential(
                r=(i + 1) * self.step_width, eps=lj_eps, sig=lj_sig, cutoff=lj_cut,
                offset=lj_off, b1=lj_b1, b2=lj_b2, e1=lj_e1, e2=lj_e2, shift=lj_shift, delta=lj_delta, lam=lj_lam)

            # Calculate forces
            f0_sim = self.system.part[0].f
            f1_sim = self.system.part[1].f
            f1_ref = self.axis * tests_common.lj_generic_force(espressomd, 
                r=(i + 1) * self.step_width, eps=lj_eps, sig=lj_sig,
                cutoff=lj_cut, offset=lj_off, b1=lj_b1, b2=lj_b2, e1=lj_e1, e2=lj_e2, delta=lj_delta, lam=lj_lam)

            # Check that energies match, ...
            self.assertFractionAlmostEqual(E_sim, E_ref)
            # force equals minus the counter-force  ...
            self.assertTrue((f0_sim == -f1_sim).all())
            # and has correct value.
            self.assertItemsFractionAlmostEqual(f1_sim, f1_ref)

        self.system.non_bonded_inter[0, 0].generic_lennard_jones.set_params(
            epsilon=0.)

    # Test Lennard-Jones Potential
    @ut.skipIf(not espressomd.has_features("LENNARD_JONES"),
               "Features not available, skipping test!")
    def test_lj(self):

        lj_eps = 1.92
        lj_sig = 1.03
        lj_cut = 1.123
        lj_shift = 0.92

        self.system.non_bonded_inter[0, 0].lennard_jones.set_params(
            epsilon=lj_eps, sigma=lj_sig, cutoff=lj_cut, shift=lj_shift)

        for i in range(113):
            self.system.part[1].pos = self.system.part[1].pos + self.step
            self.system.integrator.run(recalc_forces=True, steps=0)

            # Calculate energies
            E_sim = self.system.analysis.energy()["non_bonded"]
            E_ref = tests_common.lj_potential(
                (i + 1) * self.step_width, lj_eps, lj_sig, lj_cut, shift=lj_shift)

            # Calculate forces
            f0_sim = self.system.part[0].f
            f1_sim = self.system.part[1].f
            f1_ref = self.axis * \
                tests_common.lj_force(espressomd, r=(i + 1) * self.step_width,
                              eps=lj_eps, sig=lj_sig, cutoff=lj_cut)

            # Check that energies match, ...
            self.assertFractionAlmostEqual(E_sim, E_ref)
            # force equals minus the counter-force  ...
            self.assertTrue((f0_sim == -f1_sim).all())
            # and has correct value.
            self.assertItemsFractionAlmostEqual(f1_sim, f1_ref)

        self.system.non_bonded_inter[0, 0].lennard_jones.set_params(epsilon=0.)

    # Test Lennard-Jones Cosine Potential
    @ut.skipIf(not espressomd.has_features("LJCOS"),
               "Features not available, skipping test!")
    def test_lj_cos(self):

        ljcos_eps = 3.32
        ljcos_sig = 0.73
        ljcos_cut = 1.523
        ljcos_offset = 0.223

        self.system.non_bonded_inter[0, 0].lennard_jones_cos.set_params(
            epsilon=ljcos_eps, sigma=ljcos_sig, cutoff=ljcos_cut, offset=ljcos_offset)

        for i in range(175):
            self.system.part[1].pos = self.system.part[1].pos + self.step
            self.system.integrator.run(recalc_forces=True, steps=0)

            # Calculate energies
            E_sim = self.system.analysis.energy()["non_bonded"]
            E_ref = tests_common.lj_cos_potential(
                (i + 1) * self.step_width, eps=ljcos_eps, sig=ljcos_sig, cutoff=ljcos_cut, offset=ljcos_offset)

            # Calculate forces
            f0_sim = self.system.part[0].f
            f1_sim = self.system.part[1].f
            f1_ref = self.axis * tests_common.lj_cos_force(espressomd, (i + 1) * self.step_width,
                                                   eps=ljcos_eps, sig=ljcos_sig, cutoff=ljcos_cut, offset=ljcos_offset)

            # Check that energies match, ...
            self.assertFractionAlmostEqual(E_sim, E_ref)
            # force equals minus the counter-force  ...
            self.assertTrue((f0_sim == -f1_sim).all())
            # and has correct value.
            self.assertItemsFractionAlmostEqual(f1_sim, f1_ref)

        self.system.non_bonded_inter[
            0, 0].lennard_jones_cos.set_params(epsilon=0.)

    # Test Lennard-Jones Cosine^2 Potential
    @ut.skipIf(not espressomd.has_features("LJCOS2"),
               "Features not available, skipping test!")
    def test_lj_cos2(self):

        ljcos2_eps = 0.31
        ljcos2_sig = 0.73
        ljcos2_width = 1.523
        ljcos2_offset = 0.321

        self.system.non_bonded_inter[0, 0].lennard_jones_cos2.set_params(
            epsilon=ljcos2_eps, sigma=ljcos2_sig, offset=ljcos2_offset, width=ljcos2_width)

        for i in range(267):
            self.system.part[1].pos = self.system.part[1].pos + self.step
            self.system.integrator.run(recalc_forces=True, steps=0)

            # Calculate energies
            E_sim = self.system.analysis.energy()["non_bonded"]
            E_ref = tests_common.lj_cos2_potential(
                (i + 1) * self.step_width, eps=ljcos2_eps, sig=ljcos2_sig, offset=ljcos2_offset, width=ljcos2_width)

            # Calculate forces
            f0_sim = self.system.part[0].f
            f1_sim = self.system.part[1].f
            f1_ref = self.axis * \
                tests_common.lj_cos2_force(espressomd, 
                    r=(i + 1) * self.step_width, eps=ljcos2_eps, sig=ljcos2_sig, offset=ljcos2_offset, width=ljcos2_width)

            # Check that energies match, ...
            self.assertFractionAlmostEqual(E_sim, E_ref)
            # force equals minus the counter-force  ...
            self.assertTrue((f0_sim == -f1_sim).all())
            # and has correct value.
            self.assertItemsFractionAlmostEqual(f1_sim, f1_ref)

        self.system.non_bonded_inter[
            0, 0].lennard_jones_cos2.set_params(epsilon=0.)

    # Test Smooth-step Potential
    @ut.skipIf(not espressomd.has_features("SMOOTH_STEP"),
               "Features not available, skipping test!")
    def test_smooth_step(self):

        sst_eps = 4.92
        sst_sig = 3.03
        sst_cut = 1.253
        sst_d = 2.52
        sst_n = 11
        sst_k0 = 2.13

        self.system.non_bonded_inter[0, 0].smooth_step.set_params(
            eps=sst_eps, sig=sst_sig, cutoff=sst_cut, d=sst_d, n=sst_n, k0=sst_k0)

        for i in range(126):
            self.system.part[1].pos = self.system.part[1].pos + self.step
            self.system.integrator.run(recalc_forces=True, steps=0)

            # Calculate energies
            E_sim = self.system.analysis.energy()["non_bonded"]
            E_ref = tests_common.smooth_step_potential(
                r=(i + 1) * self.step_width, eps=sst_eps, sig=sst_sig, cutoff=sst_cut, d=sst_d, n=sst_n, k0=sst_k0)

            # Calculate forces
            f0_sim = self.system.part[0].f
            f1_sim = self.system.part[1].f
            f1_ref = self.axis * tests_common.smooth_step_force(
                r=(i + 1) * self.step_width, eps=sst_eps, sig=sst_sig, cutoff=sst_cut, d=sst_d, n=sst_n, k0=sst_k0)

            # Check that energies match, ...
            self.assertFractionAlmostEqual(E_sim, E_ref)
            # force equals minus the counter-force  ...
            self.assertTrue((f0_sim == -f1_sim).all())
            # and has correct value.
            self.assertItemsFractionAlmostEqual(f1_sim, f1_ref)

        self.system.non_bonded_inter[0, 0].smooth_step.set_params(d=0., eps=0.)

    # Test BMHTF Potential
    @ut.skipIf(not espressomd.has_features("BMHTF_NACL"),
               "Features not available, skipping test!")
    def test_bmhtf(self):

        bmhtf_a = 3.92
        bmhtf_b = 2.43
        bmhtf_c = 1.23
        bmhtf_d = 3.33
        bmhtf_sig = 0.123
        bmhtf_cut = 1.253

        self.system.non_bonded_inter[0, 0].bmhtf.set_params(
            a=bmhtf_a, b=bmhtf_b, c=bmhtf_c, d=bmhtf_d, sig=bmhtf_sig, cutoff=bmhtf_cut)

        for i in range(126):
            self.system.part[1].pos = self.system.part[1].pos + self.step
            self.system.integrator.run(recalc_forces=True, steps=0)

            # Calculate energies
            E_sim = self.system.analysis.energy()["non_bonded"]
            E_ref = tests_common.bmhtf_potential(
                r=(i + 1) * self.step_width, a=bmhtf_a, b=bmhtf_b, c=bmhtf_c, d=bmhtf_d, sig=bmhtf_sig, cutoff=bmhtf_cut)

            # Calculate forces
            f0_sim = self.system.part[0].f
            f1_sim = self.system.part[1].f
            f1_ref = self.axis * \
                tests_common.bmhtf_force(r=(i + 1) * self.step_width, a=bmhtf_a,
                                 b=bmhtf_b, c=bmhtf_c, d=bmhtf_d, sig=bmhtf_sig, cutoff=bmhtf_cut)

            # Check that energies match, ...
            self.assertFractionAlmostEqual(E_sim, E_ref)
            # force equals minus the counter-force  ...
            self.assertTrue((f0_sim == -f1_sim).all())
            # and has correct value.
            self.assertItemsFractionAlmostEqual(f1_sim, f1_ref)

        self.system.non_bonded_inter[0, 0].bmhtf.set_params(a=0., c=0., d=0.)

    # Test Morse Potential
    @ut.skipIf(not espressomd.has_features("MORSE"),
               "Features not available, skipping test!")
    def test_morse(self):

        m_eps = 1.92
        m_alpha = 3.03
        m_cut = 1.253
        m_rmin = 0.123

        self.system.non_bonded_inter[0, 0].morse.set_params(
            eps=m_eps, alpha=m_alpha, cutoff=m_cut, rmin=m_rmin)

        for i in range(126):
            self.system.part[1].pos = self.system.part[1].pos + self.step
            self.system.integrator.run(recalc_forces=True, steps=0)

            # Calculate energies
            E_sim = self.system.analysis.energy()["non_bonded"]
            E_ref = tests_common.morse_potential(
                r=(i + 1) * self.step_width, eps=m_eps, alpha=m_alpha, cutoff=m_cut, rmin=m_rmin)

            # Calculate forces
            f0_sim = self.system.part[0].f
            f1_sim = self.system.part[1].f
            f1_ref = self.axis * \
                tests_common.morse_force(r=(i + 1) * self.step_width, eps=m_eps,
                                 alpha=m_alpha, cutoff=m_cut, rmin=m_rmin)

            # Check that energies match, ...
            self.assertFractionAlmostEqual(E_sim, E_ref)
            # force equals minus the counter-force  ...
            self.assertTrue((f0_sim == -f1_sim).all())
            # and has correct value.
            self.assertItemsFractionAlmostEqual(f1_sim, f1_ref)

        self.system.non_bonded_inter[0, 0].morse.set_params(eps=0.)

    # Test Buckingham Potential
    @ut.skipIf(not espressomd.has_features("BUCKINGHAM"),
               "Features not available, skipping test!")
    def test_buckingham(self):

        b_a = 3.71
        b_b = 2.92
        b_c = 5.32
        b_d = 4.11
        b_disc = 1.03
        b_cut = 2.253
        b_shift = 0.133
        b_f1 = 0.123
        b_f2 = 0.123

        self.system.non_bonded_inter[0, 0].buckingham.set_params(
            a=b_a, b=b_b, c=b_c, d=b_d, discont=b_disc, cutoff=b_cut, shift=b_shift)

        for i in range(226):
            self.system.part[1].pos = self.system.part[1].pos + self.step
            self.system.integrator.run(recalc_forces=True, steps=0)

            # Calculate energies
            E_sim = self.system.analysis.energy()["non_bonded"]
            E_ref = tests_common.buckingham_potential(
                r=(i + 1) * self.step_width, a=b_a, b=b_b, c=b_c, d=b_d, discont=b_disc, cutoff=b_cut, shift=b_shift)

            # Calculate forces
            f0_sim = self.system.part[0].f
            f1_sim = self.system.part[1].f
            f1_ref = self.axis * \
                tests_common.buckingham_force(
                    r=(i + 1) * self.step_width, a=b_a, b=b_b, c=b_c, d=b_d, discont=b_disc, cutoff=b_cut, shift=b_shift)

            # Check that energies match, ...
            self.assertFractionAlmostEqual(E_sim, E_ref)
            # force equals minus the counter-force  ...
            self.assertTrue((f0_sim == -f1_sim).all())
            # and has correct value.
            self.assertItemsFractionAlmostEqual(f1_sim, f1_ref)

        self.system.non_bonded_inter[
            0, 0].buckingham.set_params(a=0., c=0., d=0., shift=0.)

    # Test Soft-sphere Potential
    @ut.skipIf(not espressomd.has_features("SOFT_SPHERE"),
               "Features not available, skipping test!")
    def test_soft_sphere(self):

        ss_a = 1.92
        ss_n = 3.03
        ss_cut = 1.123
        ss_off = 0.123

        self.system.non_bonded_inter[0, 0].soft_sphere.set_params(
            a=ss_a, n=ss_n, cutoff=ss_cut, offset=ss_off)

        for i in range(12):
            self.system.part[1].pos = self.system.part[1].pos + self.step
        for i in range(113):
            self.system.part[1].pos = self.system.part[1].pos + self.step
            self.system.integrator.run(recalc_forces=True, steps=0)

            # Calculate energies
            E_sim = self.system.analysis.energy()["non_bonded"]
            E_ref = tests_common.soft_sphere_potential(
                r=(i + 13) * self.step_width, a=ss_a, n=ss_n, cutoff=ss_cut, offset=ss_off)

            # Calculate forces
            f0_sim = self.system.part[0].f
            f1_sim = self.system.part[1].f
            f1_ref = self.axis * \
                tests_common.soft_sphere_force(
                    r=(i + 13) * self.step_width, a=ss_a, n=ss_n, cutoff=ss_cut, offset=ss_off)

            # Check that energies match, ...
            self.assertFractionAlmostEqual(E_sim, E_ref)
            # force equals minus the counter-force  ...
            self.assertTrue((f0_sim == -f1_sim).all())
            # and has correct value.
            self.assertItemsFractionAlmostEqual(f1_sim, f1_ref)

        self.system.non_bonded_inter[0, 0].soft_sphere.set_params(a=0.)

    # Test Hertzian Potential
    @ut.skipIf(not espressomd.has_features("HERTZIAN"),
               "Features not available, skipping test!")
    def test_hertzian(self):

        h_eps = 6.92
        h_sig = 2.432

        self.system.non_bonded_inter[0, 0].hertzian.set_params(
            eps=h_eps, sig=h_sig)

        for i in range(244):
            self.system.part[1].pos = self.system.part[1].pos + self.step
            self.system.integrator.run(recalc_forces=True, steps=0)

            # Calculate energies
            E_sim = self.system.analysis.energy()["non_bonded"]
            E_ref = tests_common.hertzian_potential(
                r=(i + 1) * self.step_width, eps=h_eps, sig=h_sig)

            # Calculate forces
            f0_sim = self.system.part[0].f
            f1_sim = self.system.part[1].f
            f1_ref = self.axis * \
                tests_common.hertzian_force(
                    r=(i + 1) * self.step_width, eps=h_eps, sig=h_sig)

            # Check that energies match, ...
            self.assertFractionAlmostEqual(E_sim, E_ref)
            # force equals minus the counter-force  ...
            self.assertTrue((f0_sim == -f1_sim).all())
            # and has correct value.
            self.assertItemsFractionAlmostEqual(f1_sim, f1_ref)

        self.system.non_bonded_inter[0, 0].hertzian.set_params(eps=0.)

    # Test Gaussian Potential
    @ut.skipIf(not espressomd.has_features("GAUSSIAN"),
               "Features not available, skipping test!")
    def test_gaussian(self):

        g_eps = 6.92
        g_sig = 4.03
        g_cut = 1.243

        self.system.non_bonded_inter[0, 0].gaussian.set_params(
            eps=g_eps, sig=g_sig, cutoff=g_cut)

        for i in range(125):
            self.system.part[1].pos = self.system.part[1].pos + self.step
            self.system.integrator.run(recalc_forces=True, steps=0)

            # Calculate energies
            E_sim = self.system.analysis.energy()["non_bonded"]
            E_ref = tests_common.gaussian_potential(
                r=(i + 1) * self.step_width, eps=g_eps, sig=g_sig, cutoff=g_cut)

            # Calculate forces
            f0_sim = self.system.part[0].f
            f1_sim = self.system.part[1].f
            f1_ref = self.axis * \
                tests_common.gaussian_force(r=(i + 1) * self.step_width,
                                    eps=g_eps, sig=g_sig, cutoff=g_cut)

            # Check that energies match, ...
            self.assertFractionAlmostEqual(E_sim, E_ref)
            # force equals minus the counter-force  ...
            self.assertTrue((f0_sim == -f1_sim).all())
            # and has correct value.
            self.assertItemsFractionAlmostEqual(f1_sim, f1_ref)

        self.system.non_bonded_inter[0, 0].gaussian.set_params(eps=0.)


if __name__ == '__main__':
    print("Features: ", espressomd.features())
    ut.main()
