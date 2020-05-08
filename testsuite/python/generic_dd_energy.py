#
# Copyright (C) 2013-2020 The ESPResSo project
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
import sys
import numpy as np
import unittest as ut
import unittest_decorators as utx
from espressomd import system, minimize_energy, generic_dd


@utx.skipIfMissingFeatures(["LENNARD_JONES"])
class Generic_DD_Energy(ut.TestCase):
    """"Test case for all generic_dd grid types.
    NVE ensemble with 10 repartitionings.
    """
    box = 64.0
    s = system.System(box_l=[box, box, box])
    partdata = box * np.random.rand(10000, 3)
    #partdata = np.loadtxt(tests_common.abspath("data/new-dd-partdata"))[:1000]
    min_partdata = None

    @classmethod
    def load_particles(cls):
        cls.s.part.clear()

        if cls.min_partdata is None:
            cls.s.part.add(pos=cls.partdata)
            minimize_energy.steepest_descent(
                cls.s, f_max=10.0, max_steps=100, gamma=0.1, max_displacement=0.01)
            cls.s.galilei.kill_particle_motion()
            cls.s.galilei.kill_particle_forces()
            cls.min_partdata = np.copy(cls.s.part[:].pos)
        else:
            cls.s.part.add(pos=cls.min_partdata)

    @classmethod
    def integrate(cls, repart_func=None):
        nsteps = 100
        for _ in range(10):
            cls.s.integrator.run(steps=nsteps // 10)
            if repart_func is not None:
                repart_func()

    @classmethod
    def get_energy(cls):
        return cls.s.analysis.energy()["total"]

    @classmethod
    def calc_ground_truth(cls):
        print("{:>15s}...".format("Ground truth"), end="")
        sys.stdout.flush()
        cls.s.cell_system.set_domain_decomposition()
        cls.load_particles()
        cls.integrate()
        cls.ground_truth = cls.get_energy()
        print(" {}".format(cls.ground_truth))

    @classmethod
    def setUpClass(cls):
        cls.s.time_step = .01
        cls.s.cell_system.skin = 0.0
        cls.s.non_bonded_inter[0, 0].lennard_jones.set_params(
            epsilon=1.0, sigma=1.0, cutoff=2.5, shift="auto")
        cls.calc_ground_truth()

    def common_test_impl(self, grid_type, skin=0.0, use_verlet=False):
        print("{:>9s}/{: 3.1f}/{:1d}...".format(grid_type, skin, use_verlet), end="")
        sys.stdout.flush()
        self.s.cell_system.skin = skin
        self.load_particles()
        gen_dd = self.s.cell_system.set_generic_dd(
            grid_type, use_verlet_lists=use_verlet)
        self.integrate(repart_func=lambda: gen_dd.repart(
            gen_dd.metric("rand")))
        en = self.get_energy()
        print(" {: 12.6f} (Rel. err.: {:.2e})".format(
            en, abs((self.ground_truth - en) / self.ground_truth)))
        self.assertAlmostEqual(self.ground_truth, en,
                               delta=abs(1e-3 * self.ground_truth))

    def common_test(self, grid_type):
        self.common_test_impl(grid_type, skin=0.0, use_verlet=False)
        self.common_test_impl(grid_type, skin=1.0, use_verlet=False)
        self.common_test_impl(grid_type, skin=0.0, use_verlet=True)
        self.common_test_impl(grid_type, skin=1.0, use_verlet=True)

    def test_all(self):
        print("Testing grid types:", generic_dd.supported_grid_types())
        for gt in generic_dd.supported_grid_types():
            self.common_test(gt)


if __name__ == "__main__":
    ut.main()
