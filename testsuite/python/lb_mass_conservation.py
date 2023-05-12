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

import espressomd
import espressomd.lb

KT = 2.25
AGRID = .5
VISC = .7
DENS = 1.7
TIME_STEP = 0.01
LB_PARAMS = {'agrid': AGRID,
             'density': DENS,
             'kinematic_viscosity': VISC,
             'tau': TIME_STEP,
             'kT': KT,
             'seed': 23}


class LBMassCommon:

    """Check the lattice-Boltzmann mass conservation."""

    system = espressomd.System(box_l=[3.0, 3.0, 3.0])
    system.time_step = TIME_STEP
    system.cell_system.skin = 0.4 * AGRID

    def setUp(self):
        self.lbf = self.lb_class(**LB_PARAMS, **self.lb_params)
        self.system.actors.add(self.lbf)
        self.system.thermostat.set_lb(LB_fluid=self.lbf, seed=3, gamma=2.0)

    def tearDown(self):
        self.system.actors.clear()
        self.system.thermostat.turn_off()

    def test_mass_conservation(self):
        self.system.integrator.run(1000)
        result = np.zeros((10, 2))
        for i in range(10):
            self.system.integrator.run(10)
            diff = self.lbf[:, :, :].density - DENS
            result[i][0] = np.mean(diff)
            result[i][1] = np.std(diff, ddof=1) / np.sqrt(np.prod(diff.shape))
        np.testing.assert_allclose(result[:, 0], 0., atol=self.atol, rtol=0.)
        np.testing.assert_array_less(result[:, 1], 0.015)


@utx.skipIfMissingFeatures(["WALBERLA"])
class LBMassWalberlaDoublePrecision(LBMassCommon, ut.TestCase):

    """Test for the Walberla implementation of the LB in double-precision."""

    lb_class = espressomd.lb.LBFluidWalberla
    lb_params = {"single_precision": False}
    atol = 1e-10


@utx.skipIfMissingFeatures(["WALBERLA"])
class LBMassWalberlaSinglePrecision(LBMassCommon, ut.TestCase):

    """Test for the Walberla implementation of the LB in single-precision."""

    lb_class = espressomd.lb.LBFluidWalberla
    lb_params = {"single_precision": True}
    atol = 5e-7


if __name__ == '__main__':
    ut.main()
