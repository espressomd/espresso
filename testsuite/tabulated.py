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
class TabulatedNonBonded(ut.TestCase):
    s = espressomd.System()
    s.box_l = 3 * [10]
    s.time_step = 0.01
    s.cell_system.skin=0.4

    def test(self):
        force = np.zeros((100,))
        energy = np.zeros((100,))
        min_ = 1.
        max_ = 2.

        dx = (max_ - min_) / 99.
        for i in range(0, 100):
            force[i] =  5 + i * 2.3 * dx
            energy[i] = 5 - i * 2.3 * dx

        self.s.non_bonded_inter[0,0].tabulated.set_params(min=min_, max=max_, energy=energy, force=force)

        self.assertTrue(np.allclose(force, self.s.non_bonded_inter[0,0].tabulated.get_params()['force']))
        self.assertTrue(np.allclose(energy, self.s.non_bonded_inter[0,0].tabulated.get_params()['energy']))
        self.assertAlmostEqual(min_, self.s.non_bonded_inter[0,0].tabulated.get_params()['min'])
        self.assertAlmostEqual(max_, self.s.non_bonded_inter[0,0].tabulated.get_params()['max'])

        self.s.part.add(id=0, type=0, pos=[5., 5., 5.0])
        self.s.part.add(id=1, type=0, pos=[5., 5., 5.5])

        # Below cutoff
        self.assertTrue(np.allclose(self.s.part[:].f, 0.0))

        for z in np.linspace(0, max_ - min_, 200, endpoint=False):
            self.s.part[1].pos = [5., 5., 6. + z]
            self.s.integrator.run(0)
            self.assertTrue(np.allclose(self.s.part[0].f, [0., 0., 5. + z * 2.3] ))
            self.assertTrue(np.allclose(self.s.part[0].f, -self.s.part[1].f))
            self.assertAlmostEqual(self.s.analysis.energy()['total'], 5. - z * 2.3)

if __name__ == "__main__":
    ut.main()
