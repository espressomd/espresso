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


class InteractionsNonBondedTest(ut.TestCase):
    system = espressomd.System()

    box_l = 10.

    start_pos = numpy.random.rand(3) * box_l
    axis = numpy.random.rand(3)
    axis /= numpy.linalg.norm(axis)
    step = axis * 0.01
    step_width = numpy.linalg.norm(step)

    def setUp(self):

        self.system.box_l = [self.box_l] * 3
        self.system.cell_system.skin = 0.4
        self.system.time_step = 1.

        self.system.part.add(id=0, pos=self.start_pos, type=0)
        self.system.part.add(id=1, pos=self.start_pos, type=0)

    def tearDown(self):

        self.system.part.clear()

    # Required, since assertAlmostEqual does NOT check significant places
    def assertFractionAlmostEqual(self, a, b, places=10):
        if abs(b) < 1E-8:
            self.assertAlmostEqual(a, b)
        else:
            self.assertAlmostEqual(a / b, 1.)

    def assertItemsFractionAlmostEqual(self, a, b):
        for i, ai in enumerate(a):
            self.assertFractionAlmostEqual(ai, b[i])

    # Analytical Expression for Generic Lennard-Jones Potential ...
    def lj_generic_potential(self, r, eps, sig, cutoff, offset=0., shift=0., e1=12., e2=6., b1=4., b2=4., delta=0., lam=1.):
        V_lj = 0.
        if (r >= offset + cutoff):
            V_lj = 0.
        else:
            # LJGEN_SOFTCORE transformations
            rroff = numpy.sqrt(
                numpy.power(r - offset, 2) + (1 - lam) * delta * sig**2)
            V_lj = eps * lam * \
                (b1 * numpy.power(sig / rroff, e1) -
                 b2 * numpy.power(sig / rroff, e2) + shift)
        return V_lj

    # ... and resulting force
    def lj_generic_force(self, r, eps, sig, cutoff, offset=0., e1=12, e2=6, b1=4., b2=4., delta=0., lam=1.):
        f_lj = 1.
        if (r >= offset + cutoff):
            f_lj = 0.
        else:
            h = (r - offset)**2 + delta * (1. - lam) * sig**2
            f_lj = (r - offset) * eps * lam * (
                b1 * e1 * numpy.power(sig / numpy.sqrt(h), e1) - b2 * e2 * numpy.power(sig / numpy.sqrt(h), e2)) / h
            if not espressomd.has_features(["LJGEN_SOFTCORE"]):
                f_lj *= numpy.sign(r - offset)
        return f_lj

    # Test Generic Lennard-Jones Potential
    @ut.skipIf(not espressomd.has_features(["LENNARD_JONES_GENERIC"]),
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
            self.system.part[1].pos += self.step
            self.system.integrator.run(recalc_forces=True, steps=0)

            # Calculate energies
            E_sim = self.system.analysis.energy()["non_bonded"]
            E_ref = self.lj_generic_potential(r=(i + 1) * self.step_width, eps=lj_eps, sig=lj_sig,
                                              cutoff=lj_cut, offset=lj_off, b1=lj_b1, b2=lj_b2, e1=lj_e1, e2=lj_e2, shift=lj_shift)

            # Calculate forces
            f0_sim = self.system.part[0].f
            f1_sim = self.system.part[1].f
            f1_ref = self.axis * self.lj_generic_force(r=(i + 1) * self.step_width, eps=lj_eps,
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
    @ut.skipIf(not espressomd.has_features(["LJGEN_SOFTCORE"]),
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
            self.system.part[1].pos += self.step
            self.system.integrator.run(recalc_forces=True, steps=0)

            # Calculate energies
            E_sim = self.system.analysis.energy()["non_bonded"]
            E_ref = self.lj_generic_potential(r=(i + 1) * self.step_width, eps=lj_eps, sig=lj_sig, cutoff=lj_cut,
                                              offset=lj_off, b1=lj_b1, b2=lj_b2, e1=lj_e1, e2=lj_e2, shift=lj_shift, delta=lj_delta, lam=lj_lam)

            # Calculate forces
            f0_sim = self.system.part[0].f
            f1_sim = self.system.part[1].f
            f1_ref = self.axis * self.lj_generic_force(r=(i + 1) * self.step_width, eps=lj_eps, sig=lj_sig,
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
    @ut.skipIf(not espressomd.has_features(["LENNARD_JONES"]),
               "Features not available, skipping test!")
    def test_lj(self):

        lj_eps = 1.92
        lj_sig = 1.03
        lj_cut = 1.123
        lj_shift = 0.92

        self.system.non_bonded_inter[0, 0].lennard_jones.set_params(
            epsilon=lj_eps, sigma=lj_sig, cutoff=lj_cut, shift=lj_shift)

        for i in range(113):
            self.system.part[1].pos += self.step
            self.system.integrator.run(recalc_forces=True, steps=0)

            # Calculate energies
            E_sim = self.system.analysis.energy()["non_bonded"]
            E_ref = self.lj_generic_potential(
                (i + 1) * self.step_width, lj_eps, lj_sig, lj_cut, shift=lj_shift * 4.)

            # Calculate forces
            f0_sim = self.system.part[0].f
            f1_sim = self.system.part[1].f
            f1_ref = self.axis * \
                self.lj_generic_force(
                    r=(i + 1) * self.step_width, eps=lj_eps, sig=lj_sig, cutoff=lj_cut)

            # Check that energies match, ...
            self.assertFractionAlmostEqual(E_sim, E_ref)
            # force equals minus the counter-force  ...
            self.assertTrue((f0_sim == -f1_sim).all())
            # and has correct value.
            self.assertItemsFractionAlmostEqual(f1_sim, f1_ref)

        self.system.non_bonded_inter[0, 0].lennard_jones.set_params(epsilon=0.)


if __name__ == '__main__':
    print("Features: ", espressomd.features())
    ut.main()
