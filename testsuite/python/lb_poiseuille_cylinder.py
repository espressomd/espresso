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

import espressomd.math
import espressomd.lb
import espressomd.observables
import espressomd.shapes

AGRID = .5
EXT_FORCE = .1
VISC = 2.7
DENS = 1.7
TIME_STEP = 0.05
BOX_L = 8.0
EFFECTIVE_RADIUS = BOX_L / 2.0 - 1.0
LB_PARAMS = {'agrid': AGRID,
             'density': DENS,
             'viscosity': VISC,
             'tau': TIME_STEP}

OBS_PARAMS = {'n_r_bins': 6,
              'n_phi_bins': 1,
              'n_z_bins': 1,
              'min_r': 0.0,
              'min_phi': -np.pi,
              'min_z': 0.0,
              'max_r': EFFECTIVE_RADIUS,
              'max_phi': np.pi,
              'max_z': BOX_L,
              'sampling_density': 1.0}


def poiseuille_flow(r, R, ext_force_density, dyn_visc):
    """
    Analytical solution for Poiseuille flow.

    Parameters
    ----------
    r : :obj:`float`
        Distance to the center of the tube.
    R : :obj:`float`
        Radius of the tube.
    ext_force_density : :obj:`float`
        Force density on the fluid parallel to the boundaries.
    dyn_visc : :obj:`float`
        Dynamic viscosity of the fluid.

    """
    return ext_force_density * 1. / (4 * dyn_visc) * (R**2.0 - r**2.0)


class LBPoiseuilleCommon:

    """
    Check the lattice-Boltzmann pressure-driven flow in a cylindrical constraint
    by comparing to the analytical solution.
    """

    system = espressomd.System(box_l=[BOX_L] * 3)
    system.time_step = TIME_STEP
    system.cell_system.skin = 0.4 * AGRID
    params = {'axis': [0, 0, 1],
              'orientation': [1, 0, 0]}

    def tearDown(self):
        self.system.actors.clear()

    def prepare(self):
        """
        Integrate the LB fluid until steady state is reached within a certain
        accuracy.

        """
        # disable periodicity except in the flow direction
        self.system.periodicity = np.logical_not(self.params['axis'])

        local_lb_params = LB_PARAMS.copy()
        local_lb_params['ext_force_density'] = np.array(
            self.params['axis']) * EXT_FORCE
        self.lbf = self.lb_class(**local_lb_params, **self.lb_params)
        self.system.actors.add(self.lbf)

        cylinder_shape = espressomd.shapes.Cylinder(
            center=self.system.box_l / 2.0, axis=self.params['axis'],
            direction=-1, radius=EFFECTIVE_RADIUS, length=BOX_L * 1.5)
        self.lbf.add_boundary_from_shape(cylinder_shape)

        # simulate until profile converges
        mid_indices = 3 * [int((BOX_L / AGRID) / 2)]
        diff = float("inf")
        old_val = self.lbf[mid_indices].velocity[2]
        while diff > 1E-5:
            self.system.integrator.run(5)
            new_val = self.lbf[mid_indices].velocity[
                np.nonzero(self.params['axis'])[0]]
            diff = abs(new_val - old_val)
            old_val = new_val

    def compare_to_analytical(self):
        """
        Compare against analytical function by calculating the RMSD.

        """
        self.prepare()
        velocities = np.zeros(int(BOX_L / AGRID))
        positions = np.zeros_like(velocities)

        for y in range(velocities.shape[0]):
            v_tmp = []
            for z in range(int(BOX_L / AGRID)):
                index = np.roll([int(BOX_L / AGRID / 2), y, z],
                                np.nonzero(self.params['axis'])[0] + 1)
                v_tmp.append(
                    self.lbf[index].velocity[np.nonzero(self.params['axis'])[0]])
            velocities[y] = np.mean(np.array(v_tmp))
            positions[y] = (y + 0.5) * AGRID

        v_measured = velocities[1:-1]
        v_expected = poiseuille_flow(
            positions[1:-1] - BOX_L / 2.0,
            EFFECTIVE_RADIUS,
            EXT_FORCE,
            VISC * DENS)
        f_half_correction = 0.5 * self.system.time_step * EXT_FORCE * AGRID**3 / DENS
        np.testing.assert_allclose(v_measured + f_half_correction,
                                   v_expected, atol=0.01, rtol=0.)

    def check_observable(self):
        if self.params['axis'] == [1, 0, 0]:
            obs_center = [0.0, BOX_L / 2.0, BOX_L / 2.0]
        elif self.params['axis'] == [0, 1, 0]:
            obs_center = [BOX_L / 2.0, 0.0, BOX_L / 2.0]
        else:
            obs_center = [BOX_L / 2.0, BOX_L / 2.0, 0.0]
        ctp = espressomd.math.CylindricalTransformationParameters(
            center=obs_center,
            axis=self.params['axis'],
            orientation=self.params['orientation'])
        local_obs_params = OBS_PARAMS.copy()
        local_obs_params['transform_params'] = ctp
        obs = espressomd.observables.CylindricalLBVelocityProfile(
            **local_obs_params)

        self.system.integrator.run(5)
        r = obs.bin_centers()[:, :, :, 0].reshape(-1)
        v_expected = poiseuille_flow(
            r,
            EFFECTIVE_RADIUS,
            EXT_FORCE,
            VISC * DENS)
        v_r, v_phi, v_z = np.copy(obs.calculate()).reshape([-1, 3]).T
        # check velocity is zero for the radial and azimuthal components
        np.testing.assert_allclose(v_r, 0., atol=1e-4, rtol=0.)
        np.testing.assert_allclose(v_phi, 0., atol=1e-4, rtol=0.)
        # check velocity is correct in the axial component
        f_half_correction = 0.5 * self.system.time_step * EXT_FORCE * AGRID**3 / DENS
        np.testing.assert_allclose(v_z + f_half_correction,
                                   v_expected, atol=3.6e-3, rtol=0.)

    def test_x(self):
        self.params['axis'] = [1, 0, 0]
        self.params['orientation'] = [0, 0, -1]
        self.compare_to_analytical()
        self.check_observable()

    def test_y(self):
        self.params['axis'] = [0, 1, 0]
        self.params['orientation'] = [1, 0, 0]
        self.compare_to_analytical()
        self.check_observable()

    def test_z(self):
        self.params['axis'] = [0, 0, 1]
        self.params['orientation'] = [1, 0, 0]
        self.compare_to_analytical()
        self.check_observable()


@utx.skipIfMissingFeatures(['WALBERLA'])
class LBPoiseuilleWalberla(LBPoiseuilleCommon, ut.TestCase):

    """Test for the Walberla implementation of the LB in double-precision."""

    lb_class = espressomd.lb.LBFluidWalberla
    lb_params = {'single_precision': False}


@utx.skipIfMissingFeatures(['WALBERLA'])
class LBPoiseuilleWalberlaSinglePrecision(LBPoiseuilleCommon, ut.TestCase):

    """Test for the Walberla implementation of the LB in single-precision."""

    lb_class = espressomd.lb.LBFluidWalberla
    lb_params = {'single_precision': True}


if __name__ == '__main__':
    ut.main()
