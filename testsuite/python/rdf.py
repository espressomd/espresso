#
# Copyright (C) 2017-2019 The ESPResSo project
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
import espressomd
import espressomd.observables
import numpy as np


class RdfTest(ut.TestCase):
    system = espressomd.System(box_l=3 * [10.])

    def tearDown(self):
        self.system.part.clear()

    def bin_volumes(self, midpoints):
        bin_size = midpoints[1] - midpoints[0]
        r = midpoints - 0.5 * bin_size
        r2 = np.power(r, 2)

        # Volumes of the bins
        return 4. * np.pi * (r2 * bin_size + r *
                             bin_size**2 + bin_size**3 / 3.)

    def test_single_type(self):
        system = self.system

        n_part = 99
        dx = system.box_l[0] / float(n_part + 1)

        for i in range(n_part):
            system.part.add(
                pos=[i * dx, 0.5 * system.box_l[1], 0.5 * system.box_l[2]], type=0)

        r_bins = 50
        r_min = 0.5 * dx
        r_max = r_bins * dx
        obs = espressomd.observables.RDF(ids1=system.part.all().id, min_r=r_min,
                                         max_r=r_max, n_r_bins=r_bins)
        rdf = obs.calculate()
        r = obs.bin_centers()
        rv = self.bin_volumes(r)
        rho = n_part / system.volume()

        parts_in_bin = rdf * rv * rho

        # All but the last bin should contain 2 particles
        np.testing.assert_allclose(parts_in_bin[:-1], 2.0, rtol=1e-1)

    def test_mixed(self):
        system = self.system

        n_part = 99
        dx = system.box_l[0] / float(n_part + 1)

        for i in range(n_part):
            system.part.add(
                pos=[i * dx, 0.5 * system.box_l[1], 0.5 * system.box_l[2]], type=(i % 2))
        partcls = system.part.all()

        r_bins = 50
        r_min = 0.5 * dx
        r_max = r_bins * dx
        obs = espressomd.observables.RDF(ids1=partcls.id[0::2],
                                         ids2=partcls.id[1::2],
                                         min_r=r_min, max_r=r_max,
                                         n_r_bins=r_bins)
        rdf01 = obs.calculate()

        r = obs.bin_centers()
        rv = self.bin_volumes(r)
        rho = 0.5 * n_part / system.volume()
        parts_in_bin = rdf01 * rv * rho

        # Every even bin should contain two parts
        np.testing.assert_allclose(parts_in_bin[0::2], 2.0, rtol=1e-1)
        # Every odd bin should contain zero parts
        np.testing.assert_allclose(parts_in_bin[1::2], 0.0)

        # Check symmetry
        obs = espressomd.observables.RDF(ids1=partcls.id[1::2],
                                         ids2=partcls.id[0::2],
                                         min_r=r_min, max_r=r_max,
                                         n_r_bins=r_bins)
        rdf10 = obs.calculate()

        np.testing.assert_allclose(rdf10, rdf01)

    def test_rdf_interface(self):
        # test setters and getters
        system = self.system
        partcls = system.part.add(pos=4 * [(0, 0, 0)], type=[0, 1, 0, 1])
        pids1 = partcls.id[0::2]
        pids2 = partcls.id[1::2]
        params = {
            'ids1': pids1,
            'ids2': pids2,
            'min_r': 1,
            'max_r': 2,
            'n_r_bins': 3}
        observable = espressomd.observables.RDF(**params)
        # check pids
        np.testing.assert_array_equal(np.copy(observable.ids1), pids1)
        np.testing.assert_array_equal(np.copy(observable.ids2), pids2)
        new_pids1 = [partcls.id[0]]
        new_pids2 = [partcls.id[1]]
        observable = espressomd.observables.RDF(
            **{**params, 'ids1': new_pids1, 'ids2': new_pids2})
        np.testing.assert_array_equal(np.copy(observable.ids1), new_pids1)
        np.testing.assert_array_equal(np.copy(observable.ids2), new_pids2)
        # check bins
        self.assertEqual(observable.n_r_bins, 3)
        observable = espressomd.observables.RDF(**{**params, 'n_r_bins': 2})
        self.assertEqual(observable.n_r_bins, 2)
        obs_data = observable.calculate()
        np.testing.assert_array_equal(obs_data.shape, [2])
        # check edges lower corner
        self.assertEqual(observable.min_r, 1)
        # check edges upper corner
        self.assertEqual(observable.max_r, 2)
        # check bin centers
        obs_bin_centers = observable.bin_centers()
        np.testing.assert_array_almost_equal(obs_bin_centers, [1.25, 1.75])
        # check exceptions
        with self.assertRaises(RuntimeError):
            espressomd.observables.RDF(**{**params, 'min_r': 100.})
        with self.assertRaises(ValueError):
            espressomd.observables.RDF(**{**params, 'n_r_bins': 0})


if __name__ == "__main__":
    ut.main()
