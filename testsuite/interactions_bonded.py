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
from tests_common import *


class InteractionsBondedTest(ut.TestCase):
    system = espressomd.System(box_l=[1.0, 1.0, 1.0])
    numpy.random.seed(seed=system.seed)

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
    def assertFractionAlmostEqual(self, a, b, places=10):
        if abs(b) < 1E-8:
            self.assertAlmostEqual(a, b)
        else:
            self.assertAlmostEqual(a / b, 1.)

    def assertItemsFractionAlmostEqual(self, a, b):
        for i, ai in enumerate(a):
            self.assertFractionAlmostEqual(ai, b[i])

    # Test Harmonic Bond
    def test_harmonic(self):

        hb_k = 5
        hb_r_0 = 1.5
        hb_r_cut = 3.355

        hb = espressomd.interactions.HarmonicBond(
            k=hb_k, r_0=hb_r_0, r_cut=hb_r_cut)
        self.system.bonded_inter.add(hb)
        self.system.part[0].add_bond((hb, 1))

        for i in range(335):
            self.system.part[1].pos = self.system.part[1].pos + self.step
            self.system.integrator.run(recalc_forces=True, steps=0)

            # Calculate energies
            E_sim = self.system.analysis.energy()["bonded"]
            E_ref = harmonic_potential(
                scalar_r=(i + 1) * self.step_width, k=hb_k, r_0=hb_r_0, r_cut=hb_r_cut)

            # Calculate forces
            f0_sim = self.system.part[0].f
            f1_sim = self.system.part[1].f
            f1_ref = self.axis * \
                harmonic_force(scalar_r=(i + 1) * self.step_width,
                                    k=hb_k, r_0=hb_r_0, r_cut=hb_r_cut)

            # Check that energies match, ...
            self.assertFractionAlmostEqual(E_sim, E_ref)
            # force equals minus the counter-force  ...
            self.assertTrue((f0_sim == -f1_sim).all())
            # and has correct value.
            self.assertItemsFractionAlmostEqual(f1_sim, f1_ref)

        # Check that bond breaks when distance > r_cut
        self.system.part[1].pos = self.system.part[1].pos + self.step
        with self.assertRaisesRegexp(Exception, "Encoutered errors during integrate"):
            self.system.integrator.run(recalc_forces=True, steps=0)

    # Test Fene Bond
    def test_fene(self):

        fene_k = 23.15
        fene_d_r_max = 3.355
        fene_r_0 = 1.1

        fene = espressomd.interactions.FeneBond(
            k=fene_k, d_r_max=fene_d_r_max, r_0=fene_r_0)
        self.system.bonded_inter.add(fene)
        self.system.part[0].add_bond((fene, 1))

        for i in range(445):
            self.system.part[1].pos = self.system.part[1].pos + self.step
            self.system.integrator.run(recalc_forces=True, steps=0)

            # Calculate energies
            E_sim = self.system.analysis.energy()["bonded"]
            E_ref = fene_potential(
                scalar_r=(i + 1) * self.step_width, k=fene_k, d_r_max=fene_d_r_max, r_0=fene_r_0)

            # Calculate forces
            f0_sim = self.system.part[0].f
            f1_sim = self.system.part[1].f
            f1_ref = self.axis * \
                fene_force(scalar_r=(i + 1) * self.step_width,
                                k=fene_k, d_r_max=fene_d_r_max, r_0=fene_r_0)

            # Check that energies match, ...
            self.assertFractionAlmostEqual(E_sim, E_ref)
            # force equals minus the counter-force  ...
            self.assertTrue((f0_sim == -f1_sim).all())
            # and has correct value.
            self.assertItemsFractionAlmostEqual(f1_sim, f1_ref)

        # Check that bond breaks when distance > r_cut
        self.system.part[1].pos = self.system.part[1].pos + self.step
        with self.assertRaisesRegexp(Exception, "Encoutered errors during integrate"):
            self.system.integrator.run(recalc_forces=True, steps=0)


if __name__ == '__main__':
    print("Features: ", espressomd.features())
    ut.main()
