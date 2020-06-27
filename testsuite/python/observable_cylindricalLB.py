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
import sys
import numpy as np
import unittest as ut
import unittest_decorators as utx
import espressomd
import espressomd.observables
import espressomd.lb
import tests_common


AGRID = 1.0
VISC = 2.7
DENS = 1.7
TIME_STEP = 0.1
LB_PARAMS = {'agrid': AGRID,
             'dens': DENS,
             'visc': VISC,
             'tau': TIME_STEP,
             }


class CylindricalLBObservableCommon:

    """
    Testcase for the CylindricalLBObservable.

    """
    lbf = None
    system = espressomd.System(box_l=(10, 10, 10))
    system.time_step = 0.01
    system.cell_system.skin = 0.4
    positions = []

    params = {
        'ids': range(10),
        'center': [5.0, 5.0, 5.0],  # center of the histogram
        'axis': 'y',
        'n_r_bins': 10,  # number of bins in r
        'n_phi_bins': 2,  # -*- in phi
        'n_z_bins': 2,  # -*- in z
        'min_r': 0.0,
        'min_phi': -np.pi,
        'min_z': -5.0,
        'max_r': 5.0,
        'max_phi': np.pi,
        'max_z': 5.0,
    }

    def tearDown(self):
        self.system.part.clear()

    def swap_axis(self, arr, axis):
        if axis == 'x':
            arr = np.dot(tests_common.rotation_matrix(
                [0, 1, 0], np.pi / 2.0), arr)
        elif axis == 'y':
            arr = np.dot(tests_common.rotation_matrix(
                [1, 0, 0], -np.pi / 2.0), arr)
        return arr

    def swap_axis_inverse(self, arr, axis):
        if axis == 'x':
            arr = np.dot(tests_common.rotation_matrix(
                [0, 1, 0], -np.pi / 2.0), arr)
        elif axis == 'y':
            arr = np.dot(tests_common.rotation_matrix(
                [1, 0, 0], np.pi / 2.0), arr)
        return arr

    def pol_coords(self):
        positions = np.zeros((len(self.positions), 3))
        for i, p in enumerate(self.positions):
            tmp = p - np.array(self.params['center'])
            tmp = self.swap_axis_inverse(tmp, self.params['axis'])
            positions[i, :] = tests_common.transform_pos_from_cartesian_to_polar_coordinates(
                tmp)
        return positions

    def set_particles(self):
        self.system.part.clear()
        self.system.part.add(pos=self.positions)

    def set_fluid_velocity(self):
        del self.positions[:]
        # Choose the cartesian velocities such that each particle gets the same
        # v_r, v_phi and v_z, respectively.
        self.v_r = .75
        self.v_phi = 2.5
        self.v_z = 1.5
        node_positions = np.arange(-4.5, 5.0, 1.0)
        for i, value in enumerate(node_positions):
            position = np.array(
                [node_positions[i], node_positions[i], node_positions[i]])
            v_y = (position[0] * np.sqrt(position[0]**2 + position[1]**2) * self.v_phi +
                   position[1] * self.v_r) / np.sqrt(position[0]**2 + position[1]**2)
            v_x = (self.v_r * np.sqrt(position[0]**2 + position[1]**2) -
                   position[1] * v_y) / position[0]
            velocity = np.array([v_x, v_y, self.v_z])
            velocity = self.swap_axis(velocity, self.params['axis'])
            position = self.swap_axis(position, self.params['axis'])
            position += np.array(self.params['center'])
            self.positions.append(position)
            self.lbf[np.array(position, dtype=int)].velocity = velocity

    def normalize_with_bin_volume(self, histogram):
        bin_volume = tests_common.get_cylindrical_bin_volume(
            self.params['n_r_bins'],
            self.params['n_phi_bins'],
            self.params['n_z_bins'],
            self.params['min_r'],
            self.params['max_r'],
            self.params['min_phi'],
            self.params['max_phi'],
            self.params['min_z'],
            self.params['max_z'])
        # Normalization
        for i in range(self.params['n_r_bins']):
            histogram[i, :, :] /= bin_volume[i]
        return histogram

    def LB_fluxdensity_profile_test(self):
        self.set_fluid_velocity()
        self.set_particles()
        # Set up the Observable.
        local_params = self.params.copy()
        if self.params['axis'] == 'x':
            local_params['axis'] = [1.0, 0.0, 0.0]
        elif self.params['axis'] == 'y':
            local_params['axis'] = [0.0, 1.0, 0.0]
        else:
            local_params['axis'] = [0.0, 0.0, 1.0]
        p = espressomd.observables.CylindricalLBFluxDensityProfileAtParticlePositions(
            **local_params)
        core_hist = np.array(p.calculate()).reshape(
            self.params['n_r_bins'],
            self.params['n_phi_bins'],
            self.params['n_z_bins'],
            3)
        core_hist_v_r = core_hist[:, :, :, 0]
        core_hist_v_phi = core_hist[:, :, :, 1]
        core_hist_v_z = core_hist[:, :, :, 2]
        self.pol_positions = self.pol_coords()
        np_hist, _ = np.histogramdd(
            self.pol_positions,
            bins=(self.params['n_r_bins'],
                  self.params['n_phi_bins'],
                  self.params['n_z_bins']),
            range=[(self.params['min_r'],
                    self.params['max_r']),
                   (self.params['min_phi'],
                    self.params['max_phi']),
                   (self.params['min_z'],
                    self.params['max_z'])])
        np_hist = self.normalize_with_bin_volume(np_hist)
        np.testing.assert_array_almost_equal(np_hist * self.v_r, core_hist_v_r)
        np.testing.assert_array_almost_equal(
            np_hist * self.v_phi, core_hist_v_phi)
        np.testing.assert_array_almost_equal(np_hist * self.v_z, core_hist_v_z)
        self.assertEqual(p.n_values(), len(np_hist.flatten()) * 3)

    def LB_velocity_profile_at_particle_positions_test(self):
        self.set_fluid_velocity()
        self.set_particles()
        # Set up the Observable.
        local_params = self.params.copy()
        if self.params['axis'] == 'x':
            local_params['axis'] = [1.0, 0.0, 0.0]
        elif self.params['axis'] == 'y':
            local_params['axis'] = [0.0, 1.0, 0.0]
        else:
            local_params['axis'] = [0.0, 0.0, 1.0]
        p = espressomd.observables.CylindricalLBVelocityProfileAtParticlePositions(
            **local_params)
        core_hist = np.array(p.calculate()).reshape(
            self.params['n_r_bins'],
            self.params['n_phi_bins'],
            self.params['n_z_bins'],
            3)
        core_hist_v_r = core_hist[:, :, :, 0]
        core_hist_v_phi = core_hist[:, :, :, 1]
        core_hist_v_z = core_hist[:, :, :, 2]
        self.pol_positions = self.pol_coords()
        np_hist, _ = np.histogramdd(
            self.pol_positions,
            bins=(self.params['n_r_bins'],
                  self.params['n_phi_bins'],
                  self.params['n_z_bins']),
            range=[(self.params['min_r'],
                    self.params['max_r']),
                   (self.params['min_phi'],
                    self.params['max_phi']),
                   (self.params['min_z'],
                    self.params['max_z'])])
        for x in np.nditer(np_hist, op_flags=['readwrite']):
            if x[...] > 0.0:
                x[...] /= x[...]
        np.testing.assert_array_almost_equal(np_hist * self.v_r, core_hist_v_r)
        np.testing.assert_array_almost_equal(
            np_hist * self.v_phi, core_hist_v_phi)
        np.testing.assert_array_almost_equal(np_hist * self.v_z, core_hist_v_z)
        self.assertEqual(p.n_values(), len(np_hist.flatten()) * 3)

    def perform_tests(self):
        self.LB_fluxdensity_profile_test()
        self.LB_velocity_profile_at_particle_positions_test()

    def test_x_axis(self):
        self.params['axis'] = 'x'
        self.perform_tests()

    def test_y_axis(self):
        self.params['axis'] = 'y'
        self.perform_tests()

    def test_z_axis(self):
        self.params['axis'] = 'z'
        self.perform_tests()

    def test_cylindrical_lb_profile_interface(self):
        # test setters and getters
        params = self.params.copy()
        params['n_r_bins'] = 4
        params['n_phi_bins'] = 6
        params['n_z_bins'] = 8
        params['axis'] = [0.0, 1.0, 0.0]
        params['sampling_density'] = 2
        del params['ids']
        observable = espressomd.observables.CylindricalLBVelocityProfile(
            **params)
        # check bins
        self.assertEqual(observable.n_r_bins, params['n_r_bins'])
        self.assertEqual(observable.n_phi_bins, params['n_phi_bins'])
        self.assertEqual(observable.n_z_bins, params['n_z_bins'])
        obs_data = observable.calculate()
        self.assertEqual(len(obs_data), 4 * 6 * 8 * 3)
        observable.n_r_bins = 1
        observable.n_phi_bins = 2
        observable.n_z_bins = 3
        self.assertEqual(observable.n_r_bins, 1)
        self.assertEqual(observable.n_phi_bins, 2)
        self.assertEqual(observable.n_z_bins, 3)
        obs_data = observable.calculate()
        self.assertEqual(len(obs_data), 1 * 2 * 3 * 3)
        # check edges lower corner
        self.assertEqual(observable.min_r, params['min_r'])
        self.assertEqual(observable.min_phi, params['min_phi'])
        self.assertEqual(observable.min_z, params['min_z'])
        observable.min_r = 4
        observable.min_phi = 5
        observable.min_z = 6
        self.assertEqual(observable.min_r, 4)
        self.assertEqual(observable.min_phi, 5)
        self.assertEqual(observable.min_z, 6)
        # check edges upper corner
        self.assertEqual(observable.max_r, params['max_r'])
        self.assertEqual(observable.max_phi, params['max_phi'])
        self.assertEqual(observable.max_z, params['max_z'])
        observable.max_r = 7
        observable.max_phi = 8
        observable.max_z = 9
        self.assertEqual(observable.max_r, 7)
        self.assertEqual(observable.max_phi, 8)
        self.assertEqual(observable.max_z, 9)
        # check center
        np.testing.assert_array_equal(
            np.copy(observable.center), params['center'])
        observable.center = [3, 2, 1]
        np.testing.assert_array_equal(np.copy(observable.center), [3, 2, 1])
        # check axis
        np.testing.assert_array_equal(np.copy(observable.axis), params['axis'])
        observable.axis = [6, 5, 4]
        np.testing.assert_array_equal(np.copy(observable.axis), [6, 5, 4])
        # check sampling_density
        self.assertEqual(observable.sampling_density, 2)
        observable.sampling_density = 3
        self.assertEqual(observable.sampling_density, 3)


class CylindricalLBObservableCPU(ut.TestCase, CylindricalLBObservableCommon):

    def setUp(self):
        self.lbf = espressomd.lb.LBFluid(**LB_PARAMS)
        self.system.actors.add(self.lbf)

    def tearDown(self):
        del self.positions[:]
        self.system.actors.remove(self.lbf)
        self.system.part.clear()


@utx.skipIfMissingGPU()
class CylindricalLBObservableGPU(ut.TestCase, CylindricalLBObservableCommon):

    def setUp(self):
        self.lbf = espressomd.lb.LBFluidGPU(**LB_PARAMS)
        self.system.actors.add(self.lbf)

    def tearDown(self):
        del self.positions[:]
        self.system.actors.remove(self.lbf)
        self.system.part.clear()


if __name__ == "__main__":
    ut.main()
