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
from espressomd import thermostat
import numpy
import unittest as ut

class InteractionsNonBondedTest(ut.TestCase):
    system = espressomd.System()

    box_l = 10.
    n_steps = 150

    start_pos = [0.123, 0.453, 1.192]
    axis = numpy.random.rand(3)
    step = axis * 3 / n_steps
    step_width = numpy.linalg.norm(step)

    # Generic Lennard-Jones Potential
    def lj_generic_potential(self, r, eps, sig, cutoff, offset=0., shift=0., e1=12, e2=6, b1=1., b2=1., delta=1., lamb=1.):
        V_lj = 0.
        if ( ( r <= offset ) or ( r >= offset + cutoff ) ):
            V_lj = 0.
        else:
            # LJGEN_SOFTCORE transformations
            rroff = numpy.sqrt( (r - offset)**2 - (1 - lamb) * delta * sig**2 )
            eps *= lamb
            V_lj = eps * ( b1 * sig**e1 / rroff**e1 - b2 * sig**e2 / rroff**e2 + shift )
        return V_lj
    
    def lj_potential(sefl, r, eps, sig, cutoff):
        V_lj = 4*eps*( sig**12 / (r-cutoff)**12 - (sig**6 / (r-cutoff)**6))
        return V_lj

    def setUp(self):

        self.system.box_l = [self.box_l]*3
        self.system.cell_system.skin = 0.4
        self.system.time_step = 1.

        self.system.part.add(id=0, pos=self.start_pos, type=0)
        self.system.part.add(id=1, pos=self.start_pos, type=0)


    def tearDown(self):

        self.system.part.clear()

    # Required, since assertAlmostEqual does NOT check significant places
    def assertFractionAlmostEqual(self, a, b, places=10):
        if b == 0:
            self.assertAlmostEqual(a, 0.)
        else:
            self.assertAlmostEqual(a/b, 1.)

    # Generic Lennard-Jones Potential
    @ut.skipIf(not espressomd.has_features(["LENNARD_JONES_GENERIC"]) ,
               "Features not available, skipping test!")
    def test_lj_generic(self):

        lj_eps = 1.0
        lj_sig = 1.0
        lj_cut = 1.12246

        self.system.non_bonded_inter[0, 0].generic_lennard_jones.set_params(
                epsilon=lj_eps, sigma=lj_sig,
                cutoff=lj_cut, shift=0,
                offset=0., e1=12, e2=6, b1=1., b2=1.)#, shift="auto")

        for i in range(self.n_steps):
            self.system.part[1].pos += self.step #self.axis * self.step_width
            self.system.integrator.run(recalc_forces=True, steps=0)

            print(self.system.part[1].pos)

            E_sim = self.system.analysis.energy()["non_bonded"]
            E_ref = self.lj_generic_potential(r=(i+1)*self.step_width, eps=lj_eps, sig=lj_sig, cutoff=lj_cut)
            print(E_sim)
            print(E_ref)
            self.assertFractionAlmostEqual(E_sim, E_ref)

        self.system.non_bonded_inter[0, 0].generic_lennard_jones.set_params(epsilon=-1)


    # Lennard-Jones Potential
    @ut.skipIf(not espressomd.has_features(["LENNARD_JONES"]) ,
               "Features not available, skipping test!")
    def test_lj(self):

        lj_eps = 1.0
        lj_sig = 1.0
        lj_cut = 1.12246

        self.system.non_bonded_inter[0, 0].lennard_jones.set_params(
            epsilon=lj_eps, sigma=lj_sig,
            cutoff=lj_cut, shift=0)#, shift="auto")

        for i in range(self.n_steps):
            self.system.part[1].pos += self.step
            self.system.integrator.run(recalc_forces=True, steps=0)

            print(self.system.part[1].pos)

            E_sim = self.system.analysis.energy()["non_bonded"]
            E_ref = self.lj_generic_potential((i+1)*self.step_width, 4*lj_eps, lj_sig, lj_cut)
            print(E_sim)
            print(E_ref)
            self.assertFractionAlmostEqual(E_sim, E_ref)

if __name__ == '__main__':
    print("Features: ", espressomd.features())
    ut.main()
