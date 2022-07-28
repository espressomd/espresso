#
# Copyright (C) 2010-2022 The ESPResSo project
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

import unittest as ut
import unittest_decorators as utx
import numpy as np
import scipy.spatial
import espressomd


@utx.skipIfMissingFeatures("LENNARD_JONES")
class AnalyzeDistributions(ut.TestCase):
    system = espressomd.System(box_l=[40., 40., 40.])
    system.cell_system.set_n_square(use_verlet_lists=False)
    np.random.seed(1234)
    # insert particles in one octant of the box to avoid PBC images
    system.part.add(pos=np.outer(np.random.random(10), system.box_l / 2.))

    def calc_min_distribution(self, bins, int_flag):
        pos = self.system.part.all().pos
        dmat = scipy.spatial.distance_matrix(pos, pos)
        dmat[np.diag_indices(dmat.shape[0])] = np.inf
        hist = np.histogram(np.min(dmat, axis=1), bins=bins, density=False)[0]
        hist = hist / np.sum(hist)
        if int_flag:
            return np.cumsum(hist)
        return hist

    def test_distribution_lin(self):
        r_min = 0.0
        r_max = 100.
        r_bins = 100
        edges = np.linspace(r_min, r_max, num=r_bins + 1, endpoint=True)
        ref_bins = (edges[1:] + edges[:-1]) / 2.
        for int_flag in (0, 1):
            ref_rdf = self.calc_min_distribution(edges, int_flag)
            core_rdf = self.system.analysis.distribution(
                type_list_a=[0], type_list_b=[0], r_min=r_min, r_max=r_max,
                r_bins=r_bins, log_flag=0, int_flag=int_flag)
            np.testing.assert_allclose(core_rdf[0], ref_bins)
            np.testing.assert_allclose(core_rdf[1], ref_rdf)

    def test_distribution_log(self):
        r_min = 0.01
        r_max = 100.
        r_bins = 100
        edges = np.geomspace(r_min, r_max, num=r_bins + 1)
        for int_flag in (0, 1):
            ref_rdf = self.calc_min_distribution(edges, int_flag)
            core_rdf = self.system.analysis.distribution(
                type_list_a=[0], type_list_b=[0], r_min=r_min, r_max=r_max,
                r_bins=r_bins, log_flag=1, int_flag=int_flag)
            ref_bins = np.geomspace(
                core_rdf[0][0], core_rdf[0][-1], r_bins, endpoint=True)
            np.testing.assert_allclose(core_rdf[0], ref_bins)
            np.testing.assert_allclose(core_rdf[1], ref_rdf)


if __name__ == "__main__":
    ut.main()
