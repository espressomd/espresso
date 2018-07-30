#
# Copyright (C) 2017 The ESPResSo project
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
import unittest as ut
import espressomd
import numpy as np

class RdfTest(ut.TestCase):
    s = espressomd.System(box_l=[1.0, 1.0, 1.0])
    s.seed = s.cell_system.get_state()['n_nodes'] * [1234]

    def setUp(self):
        self.s.box_l = 3 * [10]
        self.s.part.clear()

    def bin_volumes(self, midpoints):
        bin_size = midpoints[1] - midpoints[0]
        r = midpoints - 0.5 * bin_size
        r2 = np.power(r, 2)

        # Volumes of the bins
        return 4.*3.14159*(r2*bin_size + r * bin_size**2 + bin_size**3/3.)

    def test_single_type(self):
        s = self.s

        n_part = 99
        dx = self.s.box_l[0] / float(n_part + 1)

        for i in range(n_part):
            s.part.add(id=i, pos=[i*dx, 0.5*s.box_l[1], 0.5*s.box_l[2]], type=0)

        r_bins = 50
        r_min = 0.5 * dx
        r_max = r_bins * dx
        rdf = s.analysis.rdf(rdf_type='rdf', type_list_a=[0,1],
                             r_min=r_min, r_max=r_max, r_bins=r_bins)
        rv = self.bin_volumes(rdf[0])
        rho = n_part / (s.box_l[0]**3)

        parts_in_bin = rdf[1] * rv * rho

        # All but the last bin should contain 2 particles
        self.assertTrue(np.allclose(parts_in_bin[:-1], 2.0))

    def test_mixed(self):
        s = self.s

        n_part = 99
        dx = self.s.box_l[0] / float(n_part + 1)

        for i in range(n_part):
            s.part.add(id=i, pos=[i*dx, 0.5*s.box_l[1], 0.5*s.box_l[2]], type=(i % 2))

        r_bins = 50
        r_min = 0.5 * dx
        r_max = r_bins * dx
        rdf01 = s.analysis.rdf(rdf_type='rdf', type_list_a=[0], type_list_b=[1],
                             r_min=r_min, r_max=r_max, r_bins=r_bins)
        rv = self.bin_volumes(rdf01[0])
        rho = 0.5 * n_part / (s.box_l[0]**3)

        parts_in_bin = rdf01[1] * rv * rho

        # Every odd bin should contain two parts
        self.assertTrue(np.allclose(parts_in_bin[0:-1:2], 2.0, rtol=1e-1))
        # Every odd bin should contain zero parts
        self.assertTrue(np.allclose(parts_in_bin[1:-1:2], 0.0))

        # Check symmetry
        rdf10 = s.analysis.rdf(rdf_type='rdf', type_list_a=[1], type_list_b=[0],
                             r_min=r_min, r_max=r_max, r_bins=r_bins)

        self.assertTrue(np.allclose(rdf10, rdf01))

    def test_av(self):
        s = self.s

        for i in range(200):
            s.part.add(id = i, pos=s.box_l * np.random.random(3), type=(i % 3))

        r_bins = 50
        r_min = 0.0
        r_max = 0.49 * s.box_l[0]
        rdf = s.analysis.rdf(rdf_type='rdf', type_list_a=[0,1,2],
                             r_min=r_min, r_max=r_max, r_bins=r_bins)

        for i in range(10):
            s.analysis.append()

        rdf_av = s.analysis.rdf(rdf_type='<rdf>', type_list_a=[0,1,2],
                             r_min=r_min, r_max=r_max, r_bins=r_bins)

        self.assertTrue(np.allclose(rdf[1], rdf_av[1]))

if __name__ == "__main__":
    ut.main()
