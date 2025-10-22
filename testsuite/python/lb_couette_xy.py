#
# Copyright (C) 2021-2025 The ESPResSo project
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

import espressomd.lb
import espressomd.lees_edwards

import unittest as ut
import unittest_decorators as utx
import numpy as np


LB_PARAMS = {'agrid': 0.6,
             'density': 1.,
             'kinematic_viscosity': 1. / 6.,
             'tau': 0.5}

coord_indexes = {"x": 0, "y": 1, "z": 2}


class LBCouetteFlowCommon:

    agrid = LB_PARAMS['agrid']
    system = espressomd.System(box_l=[32 * agrid] * 3)
    system.time_step = LB_PARAMS['tau']
    system.cell_system.skin = 0.1
    system.cell_system.set_n_square()

    n_nodes = np.prod(system.cell_system.node_grid)

    def analytical(self, x, t, nu, v, h, k_max):
        """
        Analytical solution with Fourier series of the Navier-Stokes equation.

        Parameters
        ----------
        x : :obj:`float`
            Height within the channel
        t : :obj:`float`
            Time since the start up of the shear flow
        nu: :obj:`float`
            Kinematic kinematic_viscosity
        v: :obj:`float`
            Shearing velocity
        h : :obj:`float`
            Distance between shear planes
        k_max : :obj:`int`
            Upper limit of sums for sinus series

        """
        u = x / h - 0.5
        for k in np.arange(1, k_max + 1):
            wave = 2 * np.pi * k / h
            u += np.exp(-nu * wave ** 2 * t) * np.sin(wave * x) / (np.pi * k)
        return v * u

    def setUp(self):
        self.system.time = 0.
        self.system.lb = None
        self.system.lees_edwards.protocol = None

    def check_profile(self, u_getter, **kwargs):
        # carefully select the domain decomposition
        assert kwargs["shear_plane_normal"] == "y"
        h = self.system.box_l[coord_indexes[kwargs["shear_plane_normal"]]]
        agrid = self.agrid
        shear_velocity = 0.05
        k_max = 100

        protocol = espressomd.lees_edwards.LinearShear(
            shear_velocity=shear_velocity, initial_pos_offset=0., time_0=0.)
        self.system.lees_edwards.set_boundary_conditions(
            protocol=protocol, **kwargs)

        lbf = self.lb_class(**LB_PARAMS, **self.lb_params)
        self.system.lb = lbf

        # warmup
        self.system.integrator.run(16)

        # sampling
        for i in range(5, 9):
            steps = (2**i - 2**(i - 1))
            self.system.integrator.run(steps)
            pos = (np.array(range(int(h / agrid))) + 1. / 2.) * agrid
            u_ref = self.analytical(pos, self.system.time - 1., lbf.kinematic_viscosity,
                                    shear_velocity, h, k_max)
            u_lbf = np.copy(u_getter(lbf).reshape([-1]))
            np.testing.assert_allclose(u_lbf, u_ref,
                                       atol=(shear_velocity / 2.) * 1e-2, rtol=0.)

    @ut.skipIf(n_nodes == 1, "test is designed to run on multiple MPI ranks")
    @ut.expectedFailure
    def test_profile_xy_divided_shear_direction(self):
        self.system.cell_system.node_grid = [self.nodes, 1, 1]
        self.check_profile(lambda lbf: lbf[5, :, 0].velocity[:, 0],
                           shear_direction="x", shear_plane_normal="y")

    @ut.skip("TODO: LB+Lees Edwards doesn't work for domain decomposition along shear plane normal direction")  # TODO
    @ut.skipIf(n_nodes == 1, "test is designed to run on multiple MPI ranks")
    def test_profile_xy_divided_normal_direction(self):
        self.system.cell_system.node_grid = [1, self.n_nodes, 1]
        self.check_profile(lambda lbf: lbf[5, :, 0].velocity[:, 0],
                           shear_direction="x", shear_plane_normal="y")

    def test_profile_xy_divided_z_direction(self):
        self.system.cell_system.node_grid = [1, 1, self.n_nodes]
        self.check_profile(lambda lbf: lbf[5, :, 0].velocity[:, 0],
                           shear_direction="x", shear_plane_normal="y")


@utx.skipIfMissingFeatures(["WALBERLA"])
class LBCouetteFlowWalberla(LBCouetteFlowCommon, ut.TestCase):

    """Test for the Walberla implementation of the LB in double-precision."""

    lb_class = espressomd.lb.LBFluidWalberla
    lb_params = {"single_precision": False}


@utx.skipIfMissingFeatures(["WALBERLA"])
class LBCouetteFlowWalberlaSinglePrecision(LBCouetteFlowCommon, ut.TestCase):

    """Test for the Walberla implementation of the LB in single-precision."""

    lb_class = espressomd.lb.LBFluidWalberla
    lb_params = {"single_precision": True}


if __name__ == '__main__':
    ut.main()
