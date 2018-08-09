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
# Tests particle property setters/getters
from __future__ import print_function
import unittest as ut
import espressomd
import numpy as np
from time import time

@ut.skipIf(not espressomd.has_features("TABULATED"),"Skipped because feature is disabled")
class TabulatedTest(ut.TestCase):
    s = espressomd.System(box_l=[1.0, 1.0, 1.0])
    s.seed = s.cell_system.get_state()['n_nodes'] * [1234]
    s.box_l = 3 * [10]
    s.time_step = 0.01
    s.cell_system.skin=0.4

    def setUp(self):
        self.force = np.zeros((100,))
        self.energy = np.zeros((100,))
        self.min_ = 1.
        self.max_ = 2.

        self.dx = (self.max_ - self.min_) / 99.
        for i in range(0, 100):
            self.force[i] =  5 + i * 2.3 * self.dx
            self.energy[i] = 5 - i * 2.3 * self.dx

        self.s.part.clear()
        self.s.part.add(id=0, type=0, pos=[5., 5., 5.0])
        self.s.part.add(id=1, type=0, pos=[5., 5., 5.5])

    def check(self):
        # Below cutoff
        np.testing.assert_allclose(np.copy(self.s.part[:].f), 0.0)

        for z in np.linspace(0, self.max_ - self.min_, 200, endpoint=False):
            self.s.part[1].pos = [5., 5., 6. + z]
            self.s.integrator.run(0)
            np.testing.assert_allclose(np.copy(self.s.part[0].f), [0., 0., 5. + z * 2.3])
            np.testing.assert_allclose(np.copy(self.s.part[0].f), -np.copy(self.s.part[1].f))
            self.assertAlmostEqual(self.s.analysis.energy()['total'], 5. - z * 2.3)

    def test_non_bonded(self):
        self.s.non_bonded_inter[0,0].tabulated.set_params(min=self.min_, max=self.max_, energy=self.energy, force=self.force)

        np.testing.assert_allclose(self.force, self.s.non_bonded_inter[0,0].tabulated.get_params()['force'])
        np.testing.assert_allclose(self.energy, self.s.non_bonded_inter[0,0].tabulated.get_params()['energy'])
        self.assertAlmostEqual(self.min_, self.s.non_bonded_inter[0,0].tabulated.get_params()['min'])
        self.assertAlmostEqual(self.max_, self.s.non_bonded_inter[0,0].tabulated.get_params()['max'])

        self.check()

        self.s.non_bonded_inter[0,0].tabulated.set_params(min=-1, max=-1, energy=[], force=[])

    def test_bonded(self):
        from espressomd.interactions import Tabulated

        tb = Tabulated(type='distance', min=self.min_, max=self.max_, energy=self.energy, force=self.force)
        self.s.bonded_inter.add(tb)

        np.testing.assert_allclose(self.force, tb.params['force'])
        np.testing.assert_allclose(self.energy, tb.params['energy'])
        self.assertAlmostEqual(self.min_, tb.params['min'])
        self.assertAlmostEqual(self.max_, tb.params['max'])

        self.s.part[0].add_bond((tb, 1))

        self.check()

        self.s.part[0].delete_bond((tb, 1))

if __name__ == "__main__":
    ut.main()
