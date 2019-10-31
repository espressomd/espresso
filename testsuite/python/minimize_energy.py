# Copyright (C) 2010-2019 The ESPResSo project
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
import unittest as ut
import unittest_decorators as utx
import numpy as np

import espressomd


@utx.skipIfMissingFeatures("LENNARD_JONES")
class test_minimize_energy(ut.TestCase):

    np.random.seed(42)
    system = espressomd.System(box_l=[10.0, 10.0, 10.0])

    test_rotation = espressomd.has_features(("ROTATION", "DIPOLES"))
    if test_rotation:
        from espressomd.constraints import HomogeneousMagneticField

    box_l = 10.0
    density = 0.6
    vol = box_l**3
    n_part = int(vol * density)

    lj_eps = 1.0
    lj_sig = 1.0
    lj_cut = 1.12246

    def setUp(self):
        self.system.box_l = 3 * [self.box_l]
        self.system.cell_system.skin = 0.4
        self.system.time_step = 0.01
        self.system.non_bonded_inter[0, 0].lennard_jones.set_params(
            epsilon=self.lj_eps, sigma=self.lj_sig,
            cutoff=self.lj_cut, shift="auto")
        if self.test_rotation:
            self.system.constraints.add(
                self.HomogeneousMagneticField(H=[-0.5, 0, 0]))

    def tearDown(self):
        self.system.part.clear()
        self.system.integrator.set_vv()

    def test_relaxation(self):
        for i in range(self.n_part):
            p = self.system.part.add(
                id=i, pos=np.random.random(3) * self.system.box_l)
            if self.test_rotation:
                p.dip = np.random.random(3)
                p.dipm = 1
                p.rotation = (1, 1, 1)

        self.assertNotAlmostEqual(
            self.system.analysis.energy()["total"], 0, places=10)

        self.system.integrator.set_steepest_descent(
            f_max=0.0, gamma=0.1, max_displacement=0.05)

        self.system.integrator.run(500)

        self.system.constraints.clear()

        # Check
        self.assertAlmostEqual(
            self.system.analysis.energy()["total"], 0, places=10)
        np.testing.assert_allclose(np.copy(self.system.part[:].f), 0.)
        if self.test_rotation:
            np.testing.assert_allclose(np.copy(self.system.part[:].dip),
                                       np.hstack((-np.ones((self.n_part, 1)), np.zeros((self.n_part, 1)), np.zeros((self.n_part, 1)))), atol=1E-9)

    def test_rescaling(self):
        self.system.part.add(pos=[5., 5., 4.9], type=0)
        self.system.part.add(pos=[5., 5., 5.1], type=0)

        self.system.non_bonded_inter[0, 0].lennard_jones.set_params(
            epsilon=self.lj_eps, sigma=self.lj_sig,
            cutoff=self.lj_cut, shift="auto")

        self.system.integrator.run(0)
        f_old = np.copy(self.system.part[:].f)

        # No-op, because gamma = 0.
        self.system.integrator.set_steepest_descent(
            f_max=0.0, gamma=0.0, max_displacement=0.001)
        self.system.integrator.run(1)

        np.testing.assert_allclose(f_old, np.copy(self.system.part[:].f))


if __name__ == "__main__":
    ut.main()
