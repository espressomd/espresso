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
from __future__ import print_function
import unittest as ut

import espressomd

import numpy as np
import unittest_decorators as utx
required_features = ["LENNARD_JONES", "MATHEVAL"]
if espressomd.has_features(required_features):
    from espressomd.interactions import HarmonicBond, AngleHarmonic, GenericDistance, GenericAngle


@utx.skipIfMissingFeatures(["LENNARD_JONES", "MATHEVAL"])
class GenericTest(ut.TestCase):

    system = espressomd.System(box_l=[10, 10, 10])
    system.time_step = 0.01
    system.cell_system.skin = 0.1

    pos0 = system.box_l / 2 - [.1, 0, 0]
    pos1 = system.box_l / 2
    pos2 = system.box_l / 2 + [.1, 0, 0]

    system.part.add(pos=pos0, type=0)
    system.part.add(pos=pos1, type=0)
    system.part.add(pos=pos2, type=0)

    def setUp(self):
        self.system.part[0].pos = self.pos0
        self.system.part[1].pos = self.pos1
        self.system.part[2].pos = self.pos2
        self.system.part[:].v = [0, 0, 0]
        self.system.part[:].f = [0, 0, 0]

    def test_non_bonded(self):
        params = {'epsilon': 0.1, 'sigma': 0.1, 'cutoff': 1.0, 'shift': 0.0}

        # Run with LJ
        self.system.non_bonded_inter[0, 0].lennard_jones.set_params(**params)

        self.system.integrator.run(10)

        ref_pos = np.copy(self.system.part[:].pos)
        ref_nrg = self.system.analysis.energy()["non_bonded"]

        # Reset particles
        self.setUp()

        # Run with LJ from expression
        self.system.non_bonded_inter[0, 0].lennard_jones.set_params(
            **self.system.non_bonded_inter[0, 0].lennard_jones.default_params(), cutoff=0.0)

        self.system.non_bonded_inter[0, 0].generic.set_params(
            cutoff=1.0,
            energy="4*{epsilon}*(({sigma}/x)**12 - ({sigma}/x)**6)".format(**params),
            force="48*{epsilon}*(({sigma}/x)**12 - .5*({sigma}/x)**6)/x".format(**params))

        self.system.integrator.run(10)

        np.testing.assert_allclose(np.copy(self.system.part[:].pos), ref_pos)
        self.assertAlmostEqual(
            self.system.analysis.energy()["non_bonded"], ref_nrg)

    def test_bonded_distance(self):
        params = {'k': 1.0, 'r_0': 1.0, 'r_cut': 1.0}

        # Run with harmonic bond
        bond = HarmonicBond(**params)
        self.system.bonded_inter.add(bond)
        self.system.part[0].add_bond((bond, 1))

        self.system.integrator.run(10)

        ref_pos = np.copy(self.system.part[:].pos)
        ref_nrg = self.system.analysis.energy()["bonded"]

        # Reset particles
        self.setUp()
        self.system.part[0].delete_bond(self.system.part[0].bonds[0])

        # Run with LJ from expression
        bond = GenericDistance(cutoff=params['r_cut'],
                               energy="{k}/2*(x-{r_0})**2".format(**params),
                               force="-{k}*(x-{r_0})".format(**params))
        self.system.bonded_inter.add(bond)
        self.system.part[0].add_bond((bond, 1))

        self.system.integrator.run(10)

        np.testing.assert_allclose(np.copy(self.system.part[:].pos), ref_pos)
        self.assertAlmostEqual(
            self.system.analysis.energy()["bonded"], ref_nrg)

    def test_bonded_angle(self):
        params = {'bend': 1.0, 'phi0': np.pi / 2}

        # Run with harmonic bond
        bond = AngleHarmonic(**params)
        self.system.bonded_inter.add(bond)
        self.system.part[1].add_bond((bond, 0, 2))

        self.system.integrator.run(10)

        ref_pos = np.copy(self.system.part[:].pos)
        ref_nrg = self.system.analysis.energy()["bonded"]

        # Reset particles
        self.setUp()
        self.system.part[1].delete_bond(self.system.part[1].bonds[0])

        # Run with LJ from expression
        bond = GenericAngle(cutoff=np.pi,
                            energy="{bend}/2*(x-{phi0})**2".format(**params),
                            force="{bend}*(x-{phi0})".format(**params))
        self.system.bonded_inter.add(bond)
        self.system.part[1].add_bond((bond, 0, 2))

        self.system.integrator.run(10)

        np.testing.assert_allclose(np.copy(self.system.part[:].pos), ref_pos)
        self.assertAlmostEqual(
            self.system.analysis.energy()["bonded"], ref_nrg)


if __name__ == "__main__":
    ut.main()
