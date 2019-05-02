# Copyright (C) 2010-2018 The ESPResSo project
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
import numpy as np

import espressomd.lb
import espressomd.lbboundaries
import espressomd.shapes

"""
Check the Lattice Boltzmann lid driven shear flow in a slab system
by comparing to the analytical solution.

"""


AGRID = 0.6
VISC = 3.2
DENS = 2.3
TIME_STEP = 0.02
# Box size will be H +2 AGRID to make room for walls.
# The number of grid cells should be devisible by four and 3 in all directions
# for testing on multiple mpi nodes.
H = 12 * AGRID 
W = 6 * AGRID  
SHEAR_VELOCITY = 0.3

LB_PARAMS = {'agrid': AGRID,
             'dens': DENS,
             'visc': VISC,
             'tau': TIME_STEP
             }


def shear_flow(x, t, nu, v, h, k_max):
    """
    Analytical solution for driven shear flow between two plates.

    Parameters
    ----------
    x : :obj:`float`
        Position from the left plane.
    t : :obj:`float`
        Time since start of the shearing.
    nu : :obj:`float`
        Kinematic viscosity.
    v : :obj:`float`
        Shear rate.
    h : :obj:`float`
        Distance between the plates.
    k_max : :obj:`int`
        Maximum considered wave number.

    Returns
    -------
    :obj:`double` : Analytical velocity

    """

    u = x / h - 0.5
    for k in np.arange(1, k_max + 1):
        u += 1.0 / (np.pi * k) * np.exp(
            -4 * np.pi ** 2 * nu * k ** 2 / h ** 2 * t) * np.sin(2 * np.pi / h * k * x)
    return v * u


class LBShearCommon(object):

    """Base class of the test that holds the test logic."""
    lbf = None
    system = espressomd.System(box_l=[H + 2. * AGRID,
                                      W,
                                      W])
    system.time_step = TIME_STEP
    system.cell_system.skin = 0.4 * AGRID

    def check_profile(self, shear_plane_normal, shear_direction):
        """
        Integrate the LB fluid and regularly compare with
        the exact solution.

        """
        self.system.lbboundaries.clear()
        self.system.actors.clear()
        self.system.box_l = np.max(
            ((W, W, W), shear_plane_normal * (H + 2 * AGRID)), 0)

        self.system.actors.add(self.lbf)

        wall_shape1 = espressomd.shapes.Wall(
            normal=shear_plane_normal, dist=AGRID)
        wall_shape2 = espressomd.shapes.Wall(
            normal=-1.0 * shear_plane_normal, dist=-(H + AGRID))
        wall1 = espressomd.lbboundaries.LBBoundary(
            shape=wall_shape1, velocity=-0.5 * SHEAR_VELOCITY * shear_direction)
        wall2 = espressomd.lbboundaries.LBBoundary(
            shape=wall_shape2, velocity=.5 * SHEAR_VELOCITY * shear_direction)

        self.system.lbboundaries.add(wall1)
        self.system.lbboundaries.add(wall2)

        t0 = self.system.time
        sample_points = int(H / AGRID - 1)

        for i in range(9):
            self.system.integrator.run(50)

            v_expected = shear_flow(
                x=(np.arange(0, sample_points) + .5) * AGRID,
                                t=self.system.time - t0,
                                nu=VISC,
                                v=SHEAR_VELOCITY,
                                h=H,
                                k_max=100)
            for j in range(1, sample_points):
                ind = np.max(((1, 1, 1), shear_plane_normal * j + 1), 0)
                ind = np.array(ind, dtype=int)
                v_measured = self.lbf[ind[0], ind[1], ind[2]].velocity
                np.testing.assert_allclose(
                    np.copy(v_measured),
                    np.copy(v_expected[j]) * shear_direction, atol=3E-3)

        # Test stedy state stress tensor on a node
        p_eq = DENS * AGRID**2 / TIME_STEP**2 / 3
        p_expected = np.diag((p_eq, p_eq, p_eq))
        p_expected += -VISC * DENS * SHEAR_VELOCITY / H * (
            np.outer(shear_plane_normal, shear_direction)
           + np.outer(shear_direction, shear_plane_normal))
        for n in (2, 3, 4), (3, 4, 2), (5, 4, 3):
            node_stress = np.copy(self.lbf[n[0], n[1], n[2]].stress)
            np.testing.assert_allclose(node_stress,
                                       p_expected, atol=1E-5, rtol=5E-3)

        np.testing.assert_allclose(
            np.copy(wall1.get_force()),
            -np.copy(wall2.get_force()),
            atol=1E-4)
        np.testing.assert_allclose(np.copy(wall1.get_force()), 
                                   shear_direction * SHEAR_VELOCITY / H * W**2 * VISC, atol=2E-4)

    def test(self):
        x = np.array((1, 0, 0), dtype=float)
        y = np.array((0, 1, 0), dtype=float)
        z = np.array((0, 0, 1), dtype=float)
        self.check_profile(x, y)
        self.check_profile(x, z)
        self.check_profile(y, z)
        self.check_profile(x, -y)
        self.check_profile(x, -z)
        self.check_profile(y, -z)


@ut.skipIf(not espressomd.has_features(
    ['LB_BOUNDARIES']), "Skipping test due to missing features.")
class LBCPUShear(ut.TestCase, LBShearCommon):

    """Test for the CPU implementation of the LB."""

    def setUp(self):
        self.lbf = espressomd.lb.LBFluid(**LB_PARAMS)


@ut.skipIf(not espressomd.gpu_available() or not espressomd.has_features(
    ['LB_GPU', 'LB_BOUNDARIES_GPU']), "Skipping test due to missing features or gpu.")
class LBGPUShear(ut.TestCase, LBShearCommon):

    """Test for the GPU implementation of the LB."""

    def setUp(self):
        self.lbf = espressomd.lb.LBFluidGPU(**LB_PARAMS)


if __name__ == '__main__':
    ut.main()
