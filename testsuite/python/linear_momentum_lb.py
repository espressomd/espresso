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
import itertools

import espressomd.lb

"""
Check linear momentum calculation for lattice-Boltzmann.

"""


AGRID = .5
EXT_FORCE = .1
VISC = 2.7
DENS = 1.7
TIME_STEP = 0.1
BOX_L = 3.0

LB_PARAMS = {'agrid': AGRID,
             'dens': DENS,
             'visc': VISC,
             'tau': TIME_STEP,
             'ext_force_density': [0.1, 0.2, 0.3]}


class LinearMomentumTest:

    """Base class of the test that holds the test logic."""
    lbf = None
    system = espressomd.System(box_l=[BOX_L] * 3)
    system.time_step = TIME_STEP
    system.cell_system.skin = 0.4 * AGRID

    def prepare(self):
        """
        Setup random node velocities.

        """
        self.system.actors.clear()
        self.lbf = self.lbf(**LB_PARAMS)
        self.system.actors.add(self.lbf)
        for index in itertools.product(
                np.arange(0, int(np.floor(BOX_L / AGRID))), repeat=3):
            self.lbf[index].velocity = np.random.random(3) - 0.5

    def test(self):
        """
        Compare direct calculation of fluid momentum with analysis function.

        """
        self.prepare()
        linear_momentum = np.zeros(3)
        for index in itertools.product(
                np.arange(0, int(np.floor(BOX_L / AGRID))), repeat=3):
            linear_momentum += DENS * AGRID**3.0 * self.lbf[index].velocity
        analyze_linear_momentum = self.system.analysis.linear_momentum(True,  # particles
                                                                       True)  # LB fluid
        np.testing.assert_allclose(
            linear_momentum,
            analyze_linear_momentum,
            atol=1e-3)


@utx.skipIfMissingFeatures(['EXTERNAL_FORCES'])
class LBCPULinearMomentum(ut.TestCase, LinearMomentumTest):

    """Test for the CPU implementation of the LB."""

    def setUp(self):
        self.lbf = espressomd.lb.LBFluid


@utx.skipIfMissingGPU()
@utx.skipIfMissingFeatures(['LB_BOUNDARIES_GPU', 'EXTERNAL_FORCES'])
class LBGPULinearMomentum(ut.TestCase, LinearMomentumTest):

    """Test for the GPU implementation of the LB."""

    def setUp(self):
        self.lbf = espressomd.lb.LBFluidGPU


if __name__ == '__main__':
    ut.main()
