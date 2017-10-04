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

    #start_pos = [0.123, 0.453, 1.192]
    start_pos = numpy.random.rand(3)*box_l
    axis = numpy.random.rand(3)
    axis /= numpy.linalg.norm(axis)
    step = axis * 1.2 / n_steps
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
    
    def lj_potential(self, r, eps, sig, cutoff):
        V_lj = 4*eps*( sig**12 / (r-cutoff)**12 - (sig**6 / (r-cutoff)**6))
        return V_lj

    def lj_generic_force(self, r, eps, sig, cutoff, offset=0., e1=12, e2=6, b1=1., b2=1., delta=1., lamb=1.):
        f_lj = 0.
        if ( ( r <= offset ) or ( r >= offset + cutoff ) ):
            f_lj = 0.
        else:
            f_lj = (r - offset) * eps * lamb * ( (r - offset)**2 + delta * (lamb - 1) * sig**2 )**(-1. - e2/2.) * (-b2 * e2*sig + b1*e1* (sig/numpy.sqrt((r-offset)**2 + delta * (lamb - 1) * sig**2))**e1) * ( (r- offset)**2 + delta * (lamb -1) * sig**2 )*(e2/2)
            h = (r - offset)**2 + delta * (lamb - 1.) * sig**2
            f_lj = (r - offset) * eps * lamb * ( b1 * e1 * (sig / numpy.sqrt(h))**e1 - b2 * e2 * (sig / numpy.sqrt(h))**e2 ) / h 
        return f_lj

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
        if abs(b) < 1E-7 :
            self.assertAlmostEqual(a, b)
        else:
            self.assertAlmostEqual(a/b, 1.)
    
    def assertItemsFractionAlmostEqual(self, a, b):
        for i, ai in enumerate(a):
            self.assertFractionAlmostEqual(ai, b[i])

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
            self.system.part[1].pos += self.step
            self.system.integrator.run(recalc_forces=True, steps=0)

            # Calculate energies
            E_sim = self.system.analysis.energy()["non_bonded"]
            E_ref = self.lj_generic_potential(r=(i+1)*self.step_width, eps=lj_eps, sig=lj_sig, cutoff=lj_cut)

            # Calculate forces
            f0_sim = self.system.part[0].f
            f1_sim = self.system.part[1].f
            f1_ref = self.axis * self.lj_generic_force( r=(i+1)*self.step_width, eps=lj_eps, sig=lj_sig, cutoff=lj_cut)

            # Check that energies match, ...
            self.assertFractionAlmostEqual(E_sim, E_ref)
            # force and counter-force are equal, ...
            self.assertItemsEqual(f0_sim, -f1_sim)
            # and match the expected value.
            self.assertItemsFractionAlmostEqual(f1_sim, f1_ref)
             
            print(self.system.part[0].f)
            print(self.system.part[1].f)
            print(f1_ref)

            print(E_sim)
            print(E_ref)

        self.system.non_bonded_inter[0, 0].generic_lennard_jones.set_params(epsilon=0.)


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

            # Calculate energies
            E_sim = self.system.analysis.energy()["non_bonded"]
            E_ref = self.lj_generic_potential((i+1)*self.step_width, 4*lj_eps, lj_sig, lj_cut)

            # Calculate forces
            f0_sim = self.system.part[0].f
            f1_sim = self.system.part[1].f
            f1_ref = self.axis * self.lj_generic_force( r=(i+1)*self.step_width, eps=4*lj_eps, sig=lj_sig, cutoff=lj_cut)

            # Check that energies match, ...
            self.assertFractionAlmostEqual(E_sim, E_ref)
            # force and counter-force are equal, ...
            self.assertItemsEqual(f0_sim, -f1_sim)
            # and match the expected value.
            self.assertItemsFractionAlmostEqual(f1_sim, f1_ref)

        self.system.non_bonded_inter[0, 0].lennard_jones.set_params(epsilon=0.)

if __name__ == '__main__':
    print("Features: ", espressomd.features())
    ut.main()
