#
# Copyright (C) 2021-2022 The ESPResSo project
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
import numpy as np

import espressomd.lb


def u(x, t, nu, v, h, k_max):
    """
    Analytical solution with Fourier series of Navier-Stokes equation

    Parameters
    ----------
    x : :obj:`float`
        Height within the channel
    t : :obj:`float`
        Time since the start up of the shear flow
    nu: :obj:`float`
        Kinematic viscosity
    v: :obj:`float`
        Shearing velocity
    h : :obj:`float`
        Distance between shear planes
    k_max : :obj:`int`
        Upper limit of sums for sinus series
    """
    u = x / h - 0.5
    for k in np.arange(1, k_max + 1):
        u += 1.0 / (np.pi * k) * np.exp(-4 * np.pi ** 2 * nu * k ** 2 / h ** 2 * t) * np.sin(2 * np.pi / h * k * x)
    return v * u


TIME_STEP = 1.0
total_time = 1000

v = 0.0087
k_max = 100
eta = 0.27

RHO = 0.97
NU = eta / RHO
AGRID = 1.0

LB_PARAMS = {'agrid': AGRID,
             'dens': RHO,
             'visc': NU,
             'tau': TIME_STEP}


class LBCouetteFlow:

    """Base class of the test that holds the test logic."""
    h = 10.

    lbf = None
    system = espressomd.System(box_l=[h, h, h])
    system.time_step = TIME_STEP
    system.cell_system.skin = 0.4 * AGRID
    system.actors.clear()
    system.actors.add(lbf)

    for t in range(total_time):

        # Compute analytical solution
        v_expected = u(X, system.time, nu, v / time_step, box_l, k_max)

        # Read data from nodes
        X = np.arange(0, box_l) + 0.5

        # TODO finish test case for comparisson
        v_measured =

        p.testing.assert_allclose(v_measured, v_expected, atol=1e-5)

        system.integrator.run(1)


@utx.skipIfMissingFeatures("LB_WALBERLA")
class LBWalberlaPoiseuille(ut.TestCase, LBPoiseuilleCommon):

    """Test for the Walberla implementation of the LB."""

    def setUp(self):
        self.lbf = espressomd.lb.LBFluidWalberla(**LB_PARAMS)


if __name__ == '__main__':
    ut.main()
