# Copyright (C) 2010-2021 The ESPResSo project
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

import espressomd.math
import espressomd.lb
import espressomd.lbboundaries
import espressomd.observables
import espressomd.shapes
import espressomd.accumulators

AGRID = .5
VISC = 2.7
DENS = 1.7
TIME_STEP = 0.1
BOX_L = 16.0
EFFECTIVE_RADIUS = BOX_L / 2.0 - 1.0
LB_PARAMS = {'agrid': AGRID,
             'dens': DENS,
             'visc': VISC,
             'tau': TIME_STEP}

OBS_PARAMS = {'n_r_bins': 12,
              'n_phi_bins': 1,
              'n_z_bins': 1,
              'min_r': 1.0,
              'min_phi': -np.pi,
              'min_z': 0.0,
              'max_r': EFFECTIVE_RADIUS,
              'max_phi': np.pi,
              'max_z': BOX_L / 2.,
              'sampling_density': 1.0}


def taylor_couette(v1, v2, r1, r2):
    # Taylor-Couette equation
    omega1 = v1 / r1
    omega2 = v2 / r2
    eta = r1 / r2
    a = (omega2 - omega1 * eta**2) / (1. - eta**2)
    b = r1**2 * (omega1 - omega2) / (1. - eta**2)
    return a, b


class LBCircularCouetteCommon:

    """
    Check the lattice-Boltzmann velocity-driven flow in a cylindrical
    constraint by comparing to the analytical solution.
    """

    system = espressomd.System(box_l=[BOX_L, BOX_L, BOX_L / 2.])
    system.time_step = TIME_STEP
    system.cell_system.skin = 0.4 * AGRID
    params = {'axis': [0, 0, 1],
              'orientation': [1, 0, 0]}

    def tearDown(self):
        self.system.actors.clear()
        self.system.lbboundaries.clear()

    def test_taylor_couette_flow(self):
        """
        Rotate a shell filled with fluid with a non-rotating rod at the center.
        The solution to the Navier-Stokes equation, assuming an infinite rod,
        is the Taylor-Couette equation.
        """

        # disable periodicity except in the flow direction
        self.system.periodicity = np.logical_not(self.params['axis'])
        lbf = self.lb_class(**LB_PARAMS)
        self.system.actors.add(lbf)

        # create an outer cylinder that is rotating; this is achieved by
        # creating an octagon with a slip velocity parallel to each face
        sc = np.cos(np.pi / 4.)
        normals = [
            [-1, 0, 0],
            [0, -1, 0],
            [1, 0, 0],
            [0, 1, 0],
            [-sc, sc, 0],
            [sc, -sc, 0],
            [sc, sc, 0],
            [-sc, -sc, 0],
        ]
        dists = [
            2. * AGRID - BOX_L,
            2. * AGRID - BOX_L,
            2. * AGRID,
            2. * AGRID,
            2. * AGRID - BOX_L / 2.,
            2. * AGRID - BOX_L / 2.,
            2. * AGRID + BOX_L * (np.sqrt(2.) - 1.) / 2.,
            2. * AGRID - BOX_L * (1. + (np.sqrt(2.) - 1.) / 2.),
        ]
        # outer cylinder with tangential slip velocity
        slip_vel = 0.01
        for normal, dist in zip(normals, dists):
            self.system.lbboundaries.add(espressomd.lbboundaries.LBBoundary(
                shape=espressomd.shapes.Wall(normal=normal, dist=dist),
                velocity=slip_vel * np.cross(normal, self.params['axis'])))
        # inner cylinder without slip velocity
        self.system.lbboundaries.add(espressomd.lbboundaries.LBBoundary(
            shape=espressomd.shapes.Cylinder(
                center=self.system.box_l / 2.0, axis=self.params['axis'],
                direction=1, radius=1., length=BOX_L * 1.5)))

        # the system needs to be fully symmetric
        mask = np.copy(lbf[:, :, :].boundary.astype(bool))
        np.testing.assert_array_equal(mask, np.flip(mask, axis=0))
        np.testing.assert_array_equal(mask, np.flip(mask, axis=1))
        np.testing.assert_array_equal(mask, np.flip(mask, axis=2))

        # the system needs to be closed in the x and y directions
        np.testing.assert_array_equal(mask[0, :, :], 1)
        np.testing.assert_array_equal(mask[-1, :, :], 1)
        np.testing.assert_array_equal(mask[:, 0, :], 1)
        np.testing.assert_array_equal(mask[:, -1, :], 1)

        ctp = espressomd.math.CylindricalTransformationParameters(
            center=[BOX_L / 2.0, BOX_L / 2.0, 0.0],
            axis=self.params['axis'],
            orientation=self.params['orientation'])
        local_obs_params = OBS_PARAMS.copy()
        local_obs_params['transform_params'] = ctp
        obs = espressomd.observables.CylindricalLBVelocityProfile(
            **local_obs_params)

        # simulate until profile converges
        mid_indices = [int((EFFECTIVE_RADIUS / AGRID) / 2) - 2,
                       int((BOX_L / AGRID) / 2), int((BOX_L / 2. / AGRID) / 2)]
        diff = float("inf")
        old_val = lbf[mid_indices].velocity[1]
        while diff > 1e-6:
            self.system.integrator.run(10)
            new_val = lbf[mid_indices].velocity[1]
            diff = abs(new_val - old_val)
            old_val = new_val

        r = obs.bin_centers()[:, :, :, 0].reshape(-1)
        v_r, v_phi, v_z = np.copy(obs.calculate()).reshape([-1, 3]).T

        # check velocity is zero for the radial and axial components
        np.testing.assert_allclose(v_r, 0., atol=1e-6)
        np.testing.assert_allclose(v_z, 0., atol=1e-8)

        # check azimuthal velocity in the Couette regime
        a_ref, b_ref = taylor_couette(
            0.0, slip_vel, 1., BOX_L / 2. - 2. * AGRID)
        v_phi_ref = a_ref * r + b_ref / r
        v_phi_drift = np.mean(v_phi) - np.mean(v_phi_ref)
        np.testing.assert_allclose(v_phi_drift, 0., atol=5e-4)
        np.testing.assert_allclose(v_phi - v_phi_drift, v_phi_ref, atol=1e-3)


@utx.skipIfMissingFeatures(['LB_BOUNDARIES'])
class LBCPUCircularCouette(LBCircularCouetteCommon, ut.TestCase):

    """Test for the CPU implementation of the LB."""

    lb_class = espressomd.lb.LBFluid


@utx.skipIfMissingGPU()
@utx.skipIfMissingFeatures(['LB_BOUNDARIES_GPU'])
class LBGPUCircularCouette(LBCircularCouetteCommon, ut.TestCase):

    """Test for the GPU implementation of the LB."""

    lb_class = espressomd.lb.LBFluidGPU


if __name__ == '__main__':
    ut.main()
