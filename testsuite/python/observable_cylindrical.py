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
import numpy as np
import unittest as ut
import espressomd
import espressomd.observables
import espressomd.math
import tests_common


class TestCylindricalObservable(ut.TestCase):

    """
    Testcase for the cylindrical observables.

    """
    system = espressomd.System(box_l=[15.0, 15.0, 15.0])
    system.time_step = 0.01
    system.cell_system.skin = 0.4

    cyl_transform_params = espressomd.math.CylindricalTransformationParameters(
        center=3 * [7.5], axis=[1 / np.sqrt(2), 1 / np.sqrt(2), 0], orientation=[0, 0, 1])

    params = {
        'ids': None,
        'transform_params': cyl_transform_params,
        'n_r_bins': 4,
        'n_phi_bins': 3,
        'n_z_bins': 4,
        'min_r': 0.0,
        'min_phi': -np.pi,
        'min_z': -5.0,
        'max_r': 5.0,
        'max_phi': np.pi,
        'max_z': 5.0,
    }

    v_r = 0.6
    v_phi = 0.7
    v_z = 0.8

    def tearDown(self):
        self.system.part.clear()

    def calc_ellipsis_pos_vel(
            self, n_part, z_min, z_max, semi_x=1., semi_y=1.):
        """
        Calculate positions on an elliptical corkscrew line.
        Calculate cartesian velocities that lead to a
        constant velocity in cylindrical coordinates
        """

        zs = np.linspace(z_min, z_max, num=n_part)
        angles = np.linspace(-0.99 * np.pi, 0.999 * np.pi, num=n_part)

        positions = []
        velocities = []

        for angle, z in zip(angles, zs):
            position = np.array(
                [semi_x * np.cos(angle),
                 semi_y * np.sin(angle),
                 z])

            e_r, e_phi, e_z = tests_common.get_cylindrical_basis_vectors(
                position)
            velocity = self.v_r * e_r + self.v_phi * e_phi + self.v_z * e_z

            positions.append(position)
            velocities.append(velocity)

        return np.array(positions), np.array(velocities)

    def align_with_observable_frame(self, vec):
        """
        Rotate vectors from the original box frame to the frame of the observables.
        """

        # align original z to observable z
        vec = tests_common.rodrigues_rot(vec, [1, -1, 0], -np.pi / 2.)
        # original x now points along [sqrt(3),-sqrt(3),-sqrt(3)]

        # align original x to observable orientation
        vec = tests_common.rodrigues_rot(vec, [1, 1, 0], -3. / 4. * np.pi)
        return vec

    def setup_system_get_np_hist(self):
        """
        Pick positions and velocities in the original box frame
        and calculate the np histogram.
        Then rotate and move the positions and velocities
        to the frame of the observables.
        After calculating the core observables, the result should be
        the same as the np histogram obtained from the original box frame.
        """

        positions, velocities = self.calc_ellipsis_pos_vel(100, 0.99 *
                                                           self.params['min_z'], 0.9 *
                                                           self.params['max_z'], semi_x=0.9 *
                                                           self.params['max_r'], semi_y=0.2 *
                                                           self.params['max_r'])

        # first, get the numpy histogram of the cylinder coordinates
        pos_cyl = []
        for pos in positions:
            pos_cyl.append(
                tests_common.transform_pos_from_cartesian_to_polar_coordinates(pos))
        np_hist, np_edges = tests_common.get_histogram(
            np.array(pos_cyl), self.params, 'cylindrical')
        np_dens = tests_common.normalize_cylindrical_hist(
            np_hist.copy(), self.params)

        # now align the positions and velocities with the frame of reference
        # used in the observables
        pos_aligned = []
        vel_aligned = []
        for pos, vel in zip(positions, velocities):
            pos_aligned.append(
                self.align_with_observable_frame(pos) +
                self.cyl_transform_params.center)
            vel_aligned.append(self.align_with_observable_frame(vel))
        self.system.part.add(pos=pos_aligned, v=vel_aligned)
        self.params['ids'] = self.system.part.all().id

        return np_dens, np_edges

    def check_edges(self, observable, np_edges):
        core_edges = observable.call_method("edges")
        for core_edge, np_edge in zip(core_edges, np_edges):
            np.testing.assert_array_almost_equal(core_edge, np_edge)

    def test_density_profile(self):
        """
        Check that the result from the observable (in its own frame)
        matches the np result from the box frame
        """
        np_dens, np_edges = self.setup_system_get_np_hist()

        cyl_dens_prof = espressomd.observables.CylindricalDensityProfile(
            **self.params)
        core_hist = cyl_dens_prof.calculate()
        np.testing.assert_array_almost_equal(np_dens, core_hist)
        self.check_edges(cyl_dens_prof, np_edges)

    def test_vel_profile(self):
        """
        Check that the result from the observable (in its own frame)
        matches the np result from the box frame
        """
        np_dens, np_edges = self.setup_system_get_np_hist()
        cyl_vel_prof = espressomd.observables.CylindricalVelocityProfile(
            **self.params)
        core_hist = cyl_vel_prof.calculate()
        core_hist_v_r = core_hist[:, :, :, 0]
        core_hist_v_phi = core_hist[:, :, :, 1]
        core_hist_v_z = core_hist[:, :, :, 2]
        np_hist_binary = np_dens
        np_hist_binary[np.nonzero(np_hist_binary)] = 1
        np.testing.assert_array_almost_equal(
            np_hist_binary * self.v_r, core_hist_v_r)
        np.testing.assert_array_almost_equal(
            np_hist_binary * self.v_phi, core_hist_v_phi)
        np.testing.assert_array_almost_equal(
            np_hist_binary * self.v_z, core_hist_v_z)
        self.check_edges(cyl_vel_prof, np_edges)

    def test_flux_density_profile(self):
        """
        Check that the result from the observable (in its own frame)
        matches the np result from the box frame
        """
        np_dens, np_edges = self.setup_system_get_np_hist()
        cyl_flux_dens = espressomd.observables.CylindricalFluxDensityProfile(
            **self.params)
        core_hist = cyl_flux_dens.calculate()
        core_hist_v_r = core_hist[:, :, :, 0]
        core_hist_v_phi = core_hist[:, :, :, 1]
        core_hist_v_z = core_hist[:, :, :, 2]
        np.testing.assert_array_almost_equal(np_dens * self.v_r, core_hist_v_r)
        np.testing.assert_array_almost_equal(
            np_dens * self.v_phi, core_hist_v_phi)
        np.testing.assert_array_almost_equal(np_dens * self.v_z, core_hist_v_z)
        self.check_edges(cyl_flux_dens, np_edges)

    def test_cylindrical_pid_profile_interface(self):
        """
        Test setters and getters of the script interface
        """
        params = self.params.copy()
        params['n_r_bins'] = 4
        params['n_phi_bins'] = 6
        params['n_z_bins'] = 8
        self.system.part.add(pos=[0, 0, 0], type=0)
        self.system.part.add(pos=[0, 0, 0], type=1)
        params['ids'] = self.system.part.all().id
        observable = espressomd.observables.CylindricalDensityProfile(**params)
        # check pids
        np.testing.assert_array_equal(np.copy(observable.ids), params['ids'])
        with self.assertRaises(RuntimeError):
            observable.ids = [observable.ids[0]]
        # check bins
        self.assertEqual(observable.n_r_bins, params['n_r_bins'])
        self.assertEqual(observable.n_phi_bins, params['n_phi_bins'])
        self.assertEqual(observable.n_z_bins, params['n_z_bins'])
        obs_data = observable.calculate()
        np.testing.assert_array_equal(obs_data.shape, [4, 6, 8])
        observable = espressomd.observables.CylindricalDensityProfile(
            **{**params, 'n_r_bins': 1, 'n_phi_bins': 2, 'n_z_bins': 3})
        self.assertEqual(observable.n_r_bins, 1)
        self.assertEqual(observable.n_phi_bins, 2)
        self.assertEqual(observable.n_z_bins, 3)
        obs_data = observable.calculate()
        np.testing.assert_array_equal(obs_data.shape, [1, 2, 3])
        # check edges lower corner
        self.assertEqual(observable.min_r, params['min_r'])
        self.assertEqual(observable.min_phi, params['min_phi'])
        self.assertEqual(observable.min_z, params['min_z'])
        observable = espressomd.observables.CylindricalDensityProfile(
            **{**params, 'min_r': 1, 'min_phi': 2, 'min_z': 1})
        self.assertEqual(observable.min_r, 1)
        self.assertEqual(observable.min_phi, 2)
        self.assertEqual(observable.min_z, 1)
        obs_bin_edges = observable.bin_edges()
        np.testing.assert_allclose(obs_bin_edges[0, 0, 0], [1, 2, 1],
                                   rtol=0., atol=1e-14)
        # check edges upper corner
        self.assertEqual(observable.max_r, params['max_r'])
        self.assertEqual(observable.max_phi, params['max_phi'])
        self.assertEqual(observable.max_z, params['max_z'])
        observable = espressomd.observables.CylindricalDensityProfile(
            **{**params, 'max_r': 7, 'max_phi': 8, 'max_z': 9})
        self.assertEqual(observable.max_r, 7)
        self.assertEqual(observable.max_phi, 8)
        self.assertEqual(observable.max_z, 9)
        obs_bin_edges = observable.bin_edges()
        np.testing.assert_allclose(obs_bin_edges[-1, -1, -1], [7, 8, 9],
                                   rtol=0., atol=1e-14)
        # check center, axis, orientation
        ctp = espressomd.math.CylindricalTransformationParameters(
            center=[1, 2, 3], axis=[0, 1, 0], orientation=[0, 0, 1])
        observable = espressomd.observables.CylindricalDensityProfile(
            **{**params, 'transform_params': ctp})
        observable.transform_params = ctp

        for attr_name in ['center', 'axis', 'orientation']:
            np.testing.assert_array_almost_equal(np.copy(ctp.__getattr__(attr_name)),
                                                 np.copy(observable.transform_params.__getattr__(attr_name)))


if __name__ == "__main__":
    ut.main()
