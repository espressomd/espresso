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
import espressomd.observables
import tests_common


class ProfileObservablesTest(ut.TestCase):
    system = espressomd.System(box_l=[10.0, 15.0, 20.0])
    system.cell_system.skin = 0.1
    system.time_step = 0.01
    system.part.add(pos=[4.0, 4.0, 6.0], v=[0.0, 0.0, 1.0])
    system.part.add(pos=[4.0, 4.0, 6.0], v=[0.0, 0.0, 1.0])
    bin_volume = 5.0**3
    kwargs = {'ids': list(system.part.all().id),
              'n_x_bins': 2,
              'n_y_bins': 3,
              'n_z_bins': 4,
              'min_x': 0.0,
              'max_x': 10.0,
              'min_y': 0.0,
              'max_y': 15.0,
              'min_z': 0.0,
              'max_z': 20.0}

    def test_density_profile(self):
        density_profile = espressomd.observables.DensityProfile(**self.kwargs)
        obs_data = density_profile.calculate()
        obs_edges = density_profile.call_method("edges")
        obs_bin_edges = density_profile.bin_edges()
        obs_bin_centers = density_profile.bin_centers()
        np_hist, np_edges = tests_common.get_histogram(
            np.copy(self.system.part.all().pos), self.kwargs, 'cartesian',
            normed=True)
        np_hist *= len(self.system.part)
        np.testing.assert_array_almost_equal(obs_data, np_hist)
        for i in range(3):
            np.testing.assert_array_almost_equal(obs_edges[i], np_edges[i])
        np.testing.assert_array_equal(
            obs_data.shape, [self.kwargs['n_x_bins'], self.kwargs['n_y_bins'],
                             self.kwargs['n_z_bins']])
        np.testing.assert_array_almost_equal(
            obs_bin_edges[0, 0, 0],
            [self.kwargs['min_x'], self.kwargs['min_y'], self.kwargs['min_z']])
        np.testing.assert_array_almost_equal(
            obs_bin_edges[-1, -1, -1],
            [self.kwargs['max_x'], self.kwargs['max_y'], self.kwargs['max_z']])
        np.testing.assert_array_almost_equal(obs_bin_centers[0, 0, 0],
                                             [2.5, 2.5, 2.5])

    @utx.skipIfMissingFeatures("EXTERNAL_FORCES")
    def test_force_density_profile(self):
        density_profile = espressomd.observables.ForceDensityProfile(
            **self.kwargs)
        self.system.part.all().ext_force = [0.0, 0.0, 1.0]
        self.system.integrator.run(0)
        obs_data = density_profile.calculate()
        np.testing.assert_array_equal(
            obs_data.shape, [self.kwargs['n_x_bins'], self.kwargs['n_y_bins'],
                             self.kwargs['n_z_bins'], 3])
        self.assertEqual(obs_data[0, 0, 1, 2], 2.0 / self.bin_volume)
        self.assertEqual(np.sum(np.abs(obs_data)), 2.0 / self.bin_volume)

    def test_flux_density_profile(self):
        density_profile = espressomd.observables.FluxDensityProfile(
            **self.kwargs)
        self.system.integrator.run(0)
        obs_data = density_profile.calculate()
        np.testing.assert_array_equal(
            obs_data.shape, [self.kwargs['n_x_bins'], self.kwargs['n_y_bins'],
                             self.kwargs['n_z_bins'], 3])
        self.assertEqual(obs_data[0, 0, 1, 2], 2.0 / self.bin_volume)
        self.assertEqual(np.sum(np.abs(obs_data)), 2.0 / self.bin_volume)

    def test_pid_profile_interface(self):
        # test setters and getters
        params = {'ids': list(self.system.part.all().id),
                  'n_x_bins': 4,
                  'n_y_bins': 6,
                  'n_z_bins': 8,
                  'min_x': 0.0,
                  'max_x': 1.0,
                  'min_y': 2.0,
                  'max_y': 3.0,
                  'min_z': 4.0,
                  'max_z': 5.0}
        observable = espressomd.observables.DensityProfile(**params)
        # check pids
        np.testing.assert_array_equal(np.copy(observable.ids), params['ids'])
        with self.assertRaises(RuntimeError):
            observable.ids = [observable.ids[0]]
        # check bins
        self.assertEqual(observable.n_x_bins, params['n_x_bins'])
        self.assertEqual(observable.n_y_bins, params['n_y_bins'])
        self.assertEqual(observable.n_z_bins, params['n_z_bins'])
        obs_data = observable.calculate()
        np.testing.assert_array_equal(obs_data.shape, [4, 6, 8])
        observable = espressomd.observables.DensityProfile(
            **{**params, 'n_x_bins': 1, 'n_y_bins': 2, 'n_z_bins': 3})
        self.assertEqual(observable.n_x_bins, 1)
        self.assertEqual(observable.n_y_bins, 2)
        self.assertEqual(observable.n_z_bins, 3)
        obs_data = observable.calculate()
        np.testing.assert_array_equal(obs_data.shape, [1, 2, 3])
        # check edges lower corner
        self.assertEqual(observable.min_x, params['min_x'])
        self.assertEqual(observable.min_y, params['min_y'])
        self.assertEqual(observable.min_z, params['min_z'])
        observable = espressomd.observables.DensityProfile(
            **{**params, 'min_x': 0.5, 'min_y': 2, 'min_z': 1.12})
        self.assertEqual(observable.min_x, 0.5)
        self.assertEqual(observable.min_y, 2)
        self.assertEqual(observable.min_z, 1.12)
        obs_bin_edges = observable.bin_edges()
        np.testing.assert_array_almost_equal(
            obs_bin_edges[0, 0, 0], [0.5, 2.0, 1.12])
        # check edges upper corner
        self.assertEqual(observable.max_x, params['max_x'])
        self.assertEqual(observable.max_y, params['max_y'])
        self.assertEqual(observable.max_z, params['max_z'])
        observable = espressomd.observables.DensityProfile(
            **{**params, 'max_x': 7, 'max_y': 8, 'max_z': 9})
        self.assertEqual(observable.max_x, 7)
        self.assertEqual(observable.max_y, 8)
        self.assertEqual(observable.max_z, 9)
        obs_bin_edges = observable.bin_edges()
        np.testing.assert_array_equal(obs_bin_edges[-1, -1, -1], [7, 8, 9])

    def test_exceptions(self):
        params = self.kwargs.copy()
        for axis in 'xyz':
            with self.assertRaises(RuntimeError):
                espressomd.observables.DensityProfile(
                    **{**params, f'min_{axis}': 100.})
            with self.assertRaises(ValueError):
                espressomd.observables.DensityProfile(
                    **{**params, f'n_{axis}_bins': 0})


if __name__ == '__main__':
    ut.main()
