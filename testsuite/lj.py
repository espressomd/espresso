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
from tests_common import abspath

@ut.skipIf(not espressomd.has_features(["LENNARD_JONES"]) ,
           "Features not available, skipping test!")
class LennardJonesTest(ut.TestCase):
    system = espressomd.System(box_l=[1.0, 1.0, 1.0])
    data = numpy.loadtxt(abspath('data/lj_system.dat'))

    def setUp(self):
        self.system.part.clear()
        self.system.box_l = [10.7437]*3

        lj_eps = 1.0
        lj_sig = 1.0
        lj_cut = 1.12246

        self.system.non_bonded_inter[0, 0].lennard_jones.set_params(
            epsilon=lj_eps, sigma=lj_sig,
            cutoff=lj_cut, shift="auto")

        self.system.cell_system.skin = 0.4
        self.system.time_step = .1

        for i in range(self.data.shape[0]):
            self.system.part.add(id=int(self.data[i][0]), pos=[self.data[i][1], self.data[i][2], self.data[i][3]])

    def check(self):
        rms = 0.0
        max_df = 0.0

        for i in range(self.data.shape[0]):
            f = self.system.part[i].f
            for j in range(3):
                df2 = (self.data[i][4 + j] - f[j])**2
                rms += df2
                max_df = max(max_df, (df2)**0.5)

        rms = rms**0.5
        
        self.assertTrue(rms < 1e-5)
        self.assertTrue(max_df < 1e-5)

    def test_dd(self):
        self.system.cell_system.set_domain_decomposition(use_verlet_lists=False)
        self.system.integrator.run(recalc_forces=True, steps=0)

        self.check()

    def test_dd_vl(self):
        self.system.cell_system.set_domain_decomposition(use_verlet_lists=True)
        # Build VL and calc ia
        self.system.integrator.run(recalc_forces=True, steps=0)

        self.check()

        # Calc is from VLs
        self.system.integrator.run(recalc_forces=True, steps=0)
        self.check()

if __name__ == '__main__':
    print("Features: ", espressomd.features())
    ut.main()








