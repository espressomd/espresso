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


AGRID = 0.5
VISC = 1.4
DENS = 2.3
TIME_STEP = 0.01
H = 30.
SHEAR_RATE = 0.3

LB_PARAMS = {'agrid': AGRID,
             'dens': DENS,
             'visc': VISC,
             'fric': 1.0,
             'tau': TIME_STEP,
             'seed': 123}


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
                                      3.0,
                                      3.0])
    system.time_step = TIME_STEP
    system.cell_system.skin = 0.4 * AGRID

    def test_profile(self):
        """
        Integrate the LB fluid and regularly compare with
        the exact solution.

        """
        self.system.lbboundaries.clear()
        self.system.actors.clear()
        self.system.actors.add(self.lbf)

        wall_shape1 = espressomd.shapes.Wall(normal=[1, 0, 0], dist=AGRID)
        wall_shape2 = espressomd.shapes.Wall(
            normal=[-1, 0, 0], dist=-(self.system.box_l[0] - AGRID))
        wall1 = espressomd.lbboundaries.LBBoundary(
            shape=wall_shape1, velocity=[0, 0, -0.5 * SHEAR_RATE])
        wall2 = espressomd.lbboundaries.LBBoundary(
            shape=wall_shape2, velocity=[0, 0, +0.5 * SHEAR_RATE])

        self.system.lbboundaries.add(wall1)
        self.system.lbboundaries.add(wall2)

        t0 = self.system.time
        sample_points = int(H / AGRID)

        for i in range(10):
            self.system.integrator.run(100)

            v_measured = np.zeros(sample_points)
            x = np.zeros(sample_points)
            for j in range(1, sample_points + 1):
                v_measured[j - 1] = self.lbf[j, 1, 1].velocity[2]
                x[j - 1] = (j - 1 + 0.5) * AGRID

            v_expected = shear_flow(x=x,
                                    t=self.system.time - t0,
                                    nu=VISC,
                                    v=SHEAR_RATE,
                                    h=H,
                                    k_max=100)

            rmsd = np.sqrt(np.sum(np.square(v_expected - v_measured)))
            self.assertLess(rmsd, 1.e-4 * AGRID / TIME_STEP)


@ut.skipIf(not espressomd.has_features(
    ['LB', 'LB_BOUNDARIES']), "Skipping test due to missing features.")
class LBCPUShear(ut.TestCase, LBShearCommon):

    """Test for the CPU implementation of the LB."""

    def setUp(self):
        self.lbf = espressomd.lb.LBFluid(**LB_PARAMS)


@ut.skipIf(not espressomd.has_features(
    ['LB_GPU', 'LB_BOUNDARIES_GPU']), "Skipping test due to missing features.")
class LBGPUShear(ut.TestCase, LBShearCommon):

    """Test for the GPU implementation of the LB."""

    def setUp(self):
        self.lbf = espressomd.lb.LBFluidGPU(**LB_PARAMS)


if __name__ == '__main__':
    ut.main()
