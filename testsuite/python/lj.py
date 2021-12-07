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


@utx.skipIfMissingFeatures(["LENNARD_JONES"])
class LennardJonesTest(ut.TestCase):
    system = espressomd.System(box_l=[1.0, 1.0, 1.0])
    data = np.loadtxt(tests_common.abspath('data/lj_system.dat'))
    pos = data[:, 1:4]
    forces = data[:, 4:7]

    def setUp(self):
        self.system.part.clear()
        self.system.box_l = [10.7437] * 3

        lj_eps = 1.0
        lj_sig = 1.0
        lj_cut = 1.12246

        self.system.non_bonded_inter[0, 0].lennard_jones.set_params(
            epsilon=lj_eps, sigma=lj_sig, cutoff=lj_cut, shift="auto")

        self.system.cell_system.skin = 0.4
        self.system.time_step = .1

        self.system.part.add(pos=self.pos)

    def check(self):
        all_partcls = self.system.part.all()
        f_diff = np.linalg.norm(all_partcls.f - self.forces, axis=1)
        max_deviation = np.max(np.abs(all_partcls.f - self.forces))
        self.assertLess(np.mean(f_diff), 1e-7)
        self.assertLess(max_deviation, 1e-5)

    def test_dd(self):
        self.system.cell_system.set_domain_decomposition(
            use_verlet_lists=False)
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
    ut.main()
