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

import numpy as np
import unittest as ut
import unittest_decorators as utx
import espressomd
import espressomd.math
import espressomd.observables
import espressomd.lb
import tests_common


class CylindricalLBObservableCommon:

    """
    Testcase for the CylindricalLBObservables.

    """
    system = espressomd.System(box_l=3 * [14])
    system.time_step = 0.01
    system.cell_system.skin = 0.4

    lb_params = {'agrid': 1.,
                 'density': 1.2,
                 'kinematic_viscosity': 2.7,
                 'tau': 0.1,
                 }
    cyl_transform_params = espressomd.math.CylindricalTransformationParameters(
        center=3 * [7], axis=[1, 0, 0], orientation=[0, 0, 1])

    params = {
        'ids': None,
        'transform_params': cyl_transform_params,
        'n_r_bins': 4,
        'n_phi_bins': 3,
        'n_z_bins': 5,
        'min_r': 0.0,
        'min_phi': -np.pi,
        'min_z': -6.0,
        'max_r': 6.0,
        'max_phi': np.pi,
        'max_z': 6.0,
    }

    v_r = 0.02
    v_phi = 0.04
    v_z = 0.03

    def setUp(self):
        self.lbf = self.lb_class(**self.lb_params, **self.lb_params_extra)
        self.system.actors.add(self.lbf)

    def tearDown(self):
        self.system.actors.clear()
        self.system.part.clear()

    def calc_vel_at_pos(self, positions):
        """
        In cylindrical coordinates, all velocities are the same.
        In cartesian they depend on the position.
        The cartesian velocities are calculated here.
        """

        vels = []
        for pos in positions:
            e_r, e_phi, e_z = tests_common.get_cylindrical_basis_vectors(pos)
            velocity = self.v_r * e_r + self.v_phi * e_phi + self.v_z * e_z
            vels.append(velocity)
        return vels

    def align_with_observable_frame(self, vec):
        """
        Rotate vectors from the original box frame to
        the frame of the observables.
        """

        # align original z to observable z
        vec = tests_common.rodrigues_rot(vec, [0, 1, 0], np.pi / 2.)
        # original x now points along [0,0,-1]

        # align original x to observable orientation
        vec = tests_common.rodrigues_rot(vec, [1, 0, 0], np.pi)
        return vec

    def setup_system_get_np_hist(self):
        """
        Pick positions and velocities in the original box frame and
        calculate the np histogram. Then rotate and move the positions
        and velocities to the frame of the observables.
        After calculating the core observables, the result should be
        the same as the np histogram obtained from the original box frame.
        """

        nodes = np.array(np.meshgrid([1, 2], [1, 2], [
                         1, 1, 1, 1, 2])).T.reshape(-1, 3)
        positions = nodes + 3 * [0.5]
        velocities = self.calc_vel_at_pos(positions)

        # get the histogram from numpy
        pos_cyl = []
        for pos in positions:
            pos_cyl.append(
                tests_common.transform_pos_from_cartesian_to_polar_coordinates(pos))
        np_hist, np_edges = tests_common.get_histogram(
            np.array(pos_cyl), self.params, 'cylindrical')

        # the particles only determine the evaluation points, not the values of
        # the observables
        np_hist[np.nonzero(np_hist)] = 1

        # now align the positions and velocities with the frame of reference
        # used in the observables
        pos_aligned = []
        vel_aligned = []
        for pos, vel in zip(positions, velocities):
            pos_aligned.append(
                self.align_with_observable_frame(pos) +
                self.cyl_transform_params.center)
            vel_aligned.append(self.align_with_observable_frame(vel))
        node_aligned = np.array(
            np.rint(
                np.array(pos_aligned) -
                3 *
                [0.5]),
            dtype=int)
        self.system.part.add(pos=pos_aligned, v=vel_aligned)
        self.params['ids'] = self.system.part.all().id

        for node, vel in zip(node_aligned, vel_aligned):
            self.lbf[node].velocity = vel

        return np_hist, np_edges

    def check_edges(self, observable, np_edges):
        core_edges = observable.call_method("edges")
        for core_edge, np_edge in zip(core_edges, np_edges):
            np.testing.assert_array_almost_equal(core_edge, np_edge)

    def test_cylindrical_lb_vel_profile_obs(self):
        """
        Check that the result from the observable (in its own frame)
        matches the np result from the box frame
        """

        np_hist_binary, np_edges = self.setup_system_get_np_hist()
        vel_obs = espressomd.observables.CylindricalLBVelocityProfileAtParticlePositions(
            **self.params)
        core_hist_v = vel_obs.calculate()
        core_hist_v_r = core_hist_v[:, :, :, 0]
        core_hist_v_phi = core_hist_v[:, :, :, 1]
        core_hist_v_z = core_hist_v[:, :, :, 2]
        np.testing.assert_array_almost_equal(
            np_hist_binary * self.v_r, core_hist_v_r)
        np.testing.assert_array_almost_equal(
            np_hist_binary * self.v_phi, core_hist_v_phi)
        np.testing.assert_array_almost_equal(
            np_hist_binary * self.v_z, core_hist_v_z)
        self.check_edges(vel_obs, np_edges)

    def test_cylindrical_lb_profile_interface(self):
        """
        Test setters and getters of the script interface
        """

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
        np.testing.assert_array_equal(obs_data.shape, [4, 6, 8, 3])
        observable = espressomd.observables.CylindricalLBVelocityProfile(
            **{**params, 'n_r_bins': 1, 'n_phi_bins': 2, 'n_z_bins': 3})
        self.assertEqual(observable.n_r_bins, 1)
        self.assertEqual(observable.n_phi_bins, 2)
        self.assertEqual(observable.n_z_bins, 3)
        obs_data = observable.calculate()
        np.testing.assert_array_equal(obs_data.shape, [1, 2, 3, 3])
        # check edges lower corner
        self.assertEqual(observable.min_r, params['min_r'])
        self.assertEqual(observable.min_phi, params['min_phi'])
        self.assertEqual(observable.min_z, params['min_z'])
        observable = espressomd.observables.CylindricalLBVelocityProfile(
            **{**params, 'min_r': 4.0, 'min_phi': 1.0, 'min_z': 1.0})
        self.assertEqual(observable.min_r, 4.0)
        self.assertEqual(observable.min_phi, 1.0)
        self.assertEqual(observable.min_z, 1.0)
        obs_bin_edges = observable.bin_edges()
        np.testing.assert_allclose(obs_bin_edges[0, 0, 0], [4, 1, 1],
                                   rtol=0., atol=1e-14)
        # check edges upper corner
        self.assertEqual(observable.max_r, params['max_r'])
        self.assertEqual(observable.max_phi, params['max_phi'])
        self.assertEqual(observable.max_z, params['max_z'])
        observable = espressomd.observables.CylindricalLBVelocityProfile(
            **{**params, 'max_r': 7, 'max_phi': 8, 'max_z': 9})
        self.assertEqual(observable.max_r, 7)
        self.assertEqual(observable.max_phi, 8)
        self.assertEqual(observable.max_z, 9)
        obs_bin_edges = observable.bin_edges()
        np.testing.assert_allclose(obs_bin_edges[-1, -1, -1], [7, 8, 9],
                                   rtol=0., atol=1e-14)
        # check sampling_density
        self.assertEqual(
            observable.sampling_density,
            params['sampling_density'])
        observable = espressomd.observables.CylindricalLBVelocityProfile(
            **{**params, 'sampling_density': 3})
        self.assertEqual(observable.sampling_density, 3)
        # check center, axis, orientation
        ctp = espressomd.math.CylindricalTransformationParameters(
            center=[1, 2, 3], axis=[0, 1, 0], orientation=[0, 0, 1])
        observable = espressomd.observables.CylindricalLBVelocityProfile(
            **{**params, 'transform_params': ctp})
        observable.transform_params = ctp

        for attr_name in ['center', 'axis', 'orientation']:
            np.testing.assert_array_almost_equal(np.copy(ctp.__getattr__(attr_name)),
                                                 np.copy(observable.transform_params.__getattr__(attr_name)))

    def test_cylindrical_lb_flux_density_obs(self):
        """
        Check that the result from the observable (in its own frame)
        matches the np result from the box frame.
        Only for CPU because density interpolation is not implemented for GPU LB.
        """
        np_hist_binary, np_edges = self.setup_system_get_np_hist()

        flux_obs = espressomd.observables.CylindricalLBFluxDensityProfileAtParticlePositions(
            **self.params)
        core_hist_fl = flux_obs.calculate()
        core_hist_fl_r = core_hist_fl[:, :, :, 0]
        core_hist_fl_phi = core_hist_fl[:, :, :, 1]
        core_hist_fl_z = core_hist_fl[:, :, :, 2]

        np.testing.assert_array_almost_equal(
            np_hist_binary *
            self.lb_params['density'] *
            self.v_r,
            core_hist_fl_r)
        np.testing.assert_array_almost_equal(
            np_hist_binary *
            self.lb_params['density'] *
            self.v_phi,
            core_hist_fl_phi)
        np.testing.assert_array_almost_equal(
            np_hist_binary *
            self.lb_params['density'] *
            self.v_z,
            core_hist_fl_z)
        self.check_edges(flux_obs, np_edges)


@utx.skipIfMissingFeatures(["WALBERLA"])
class CylindricalLBObservableWalberla(
        CylindricalLBObservableCommon, ut.TestCase):

    """Test for the Walberla implementation of the LB in double-precision."""

    lb_class = espressomd.lb.LBFluidWalberla
    lb_params_extra = {"single_precision": False}


@utx.skipIfMissingFeatures(["WALBERLA"])
class CylindricalLBObservableWalberlaSinglePrecision(
        CylindricalLBObservableWalberla, ut.TestCase):

    """Test for the Walberla implementation of the LB in single-precision."""

    lb_class = espressomd.lb.LBFluidWalberla
    lb_params_extra = {"single_precision": True}


if __name__ == "__main__":
    ut.main()
