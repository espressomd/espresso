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
        for p in range(cls.num_part):
            cls.system.part.add(
                id=p,
                pos=np.random.random() * cls.system.box_l)

    def calc_rdf(self, r, bins):
        # this generates indices for all i<j combinations
        ij = np.triu_indices(len(r), k=1)
        r_ij = r[ij[0]] - r[ij[1]]
        dist = np.sqrt(np.sum(r_ij**2, axis=1))
        hist = np.histogram(dist, bins=bins, density=False)[0]
        return hist

    def calc_min_distribution(self, bins):
        dist = []
        for i in range(self.num_part):
            dist.append(np.min([self.system.distance(
                self.system.part[i], p.pos) for p in self.system.part if p.id != i]))
        hist = np.histogram(dist, bins=bins, density=False)[0]
        return hist / (float(np.sum(hist)))

    # test system.analysis.rdf()
    def test_rdf(self):
        # increase PBC to remove mirror images
        old_pos = self.system.part[:].pos.copy()
        self.system.box_l = self.system.box_l * 2.
        self.system.part[:].pos = old_pos
        r_min = 0.0
        r_max = 100.0
        r_bins = 10
        bin_width = (r_max - r_min) / r_bins
        bins = np.arange(r_min, r_max + bin_width, bin_width)
        bin_volume = 4. / 3. * np.pi * (bins[1:]**3 - bins[:-1]**3)
        box_volume = np.prod(np.copy(self.system.box_l))
        # all the same type
        core_rdf = self.system.analysis.rdf(rdf_type='rdf',
                                            type_list_a=[0],
                                            type_list_b=[0],
                                            r_min=r_min,
                                            r_max=r_max,
                                            r_bins=r_bins)
        num_pair = 0.5 * (self.num_part) * (self.num_part - 1)
        r = self.system.part[:].pos
        # bins
        np.testing.assert_allclose(core_rdf[0], (bins[1:] + bins[:-1]) * 0.5)
        # rdf
        np.testing.assert_allclose(
            core_rdf[1] * bin_volume * num_pair / box_volume,
            self.calc_rdf(r, bins))
        # change one type
        self.system.part[0].type = 1
        r = self.system.part[1:].pos
        core_rdf = self.system.analysis.rdf(rdf_type='rdf',
                                            type_list_a=[0],
                                            type_list_b=[0],
                                            r_min=r_min,
                                            r_max=r_max,
                                            r_bins=r_bins)
        num_pair = 0.5 * (self.num_part - 1) * (self.num_part - 2)
        np.testing.assert_allclose(
            core_rdf[1] * bin_volume * num_pair / box_volume,
            self.calc_rdf(r, bins))

        # compare with type
        core_rdf = self.system.analysis.rdf(rdf_type='rdf',
                                            type_list_a=[1],
                                            type_list_b=[0],
                                            r_min=r_min,
                                            r_max=r_max,
                                            r_bins=r_bins)
        num_pair = (self.num_part - 1)
        dist = np.sqrt(
            np.sum((self.system.part[1:].pos - self.system.part[0].pos)**2, axis=1))
        hist = np.histogram(dist, bins=bins, density=False)[0]
        np.testing.assert_allclose(
            core_rdf[1] * bin_volume * num_pair / box_volume, hist)
        # restore PBC
        self.system.box_l = self.system.box_l / 2.
        self.system.part[:].pos = old_pos

    # test system.analysis.distribution(), all the same particle types
    def test_distribution_lin(self):
        # increase PBC to remove mirror images
        old_pos = self.system.part[:].pos.copy()
        self.system.box_l = self.system.box_l * 2.
        self.system.part[:].pos = old_pos
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
