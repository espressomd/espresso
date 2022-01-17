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
class AnalyzeDistributions(ut.TestCase):
    system = espressomd.System(box_l=[1.0, 1.0, 1.0])
    np.random.seed(1234)
    num_part = 10

    @classmethod
    def setUpClass(cls):
        box_l = 20.0
        # start with a small box
        cls.system.box_l = np.array([box_l, box_l, box_l])
        cls.system.cell_system.set_n_square(use_verlet_lists=False)
        cls.partcls = cls.system.part.add(
            pos=np.outer(np.random.random(cls.num_part), cls.system.box_l))

    def calc_min_distribution(self, bins):
        dist = []
        for p1 in self.system.part:
            dist.append(min(self.system.distance(p1, p2.pos)
                            for p2 in self.system.part if p1.id != p2.id))
        hist = np.histogram(dist, bins=bins, density=False)[0]
        return hist / (float(np.sum(hist)))

    # test system.analysis.distribution(), all the same particle types
    def test_distribution_lin(self):
        # increase PBC to remove mirror images
        old_pos = self.partcls.pos.copy()
        self.system.box_l = self.system.box_l * 2.
        self.partcls.pos = old_pos
        r_min = 0.0
        r_max = 100.0
        r_bins = 100
        bins = np.linspace(r_min, r_max, num=r_bins + 1, endpoint=True)
        # no int flag
        core_rdf = self.system.analysis.distribution(type_list_a=[0],
                                                     type_list_b=[0],
                                                     r_min=r_min,
                                                     r_max=r_max,
                                                     r_bins=r_bins,
                                                     log_flag=0,
                                                     int_flag=0)
        # bins
        np.testing.assert_allclose(core_rdf[0], (bins[1:] + bins[:-1]) * 0.5)

        # rdf
        np.testing.assert_allclose(core_rdf[1],
                                   self.calc_min_distribution(bins))
        # with int flag
        core_rdf = self.system.analysis.distribution(type_list_a=[0],
                                                     type_list_b=[0],
                                                     r_min=r_min,
                                                     r_max=r_max,
                                                     r_bins=r_bins,
                                                     log_flag=0,
                                                     int_flag=1)
        np.testing.assert_allclose(core_rdf[1],
                                   np.cumsum(self.calc_min_distribution(bins)))


if __name__ == "__main__":
    ut.main()
