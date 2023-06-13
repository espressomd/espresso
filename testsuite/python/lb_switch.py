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
import itertools


@utx.skipIfMissingFeatures(["EXTERNAL_FORCES"])
class LBSwitchActor(ut.TestCase):
    system = espressomd.System(box_l=[10.0, 10.0, 10.0])

    system.time_step = 0.01
    system.cell_system.skin = 0.1

    def switch_test(self, GPU=False):
        system = self.system
        system.actors.clear()
        p = system.part.add(pos=[1., 1., 1.], v=[1., 0, 0], fix=3 * [True])

        lb_fluid_params = {'agrid': 2.0, 'dens': 1.0, 'visc': 1.0, 'tau': 0.03}
        friction_1 = 1.5
        friction_2 = 4.0

        if GPU:
            lb_fluid_1 = espressomd.lb.LBFluidGPU(**lb_fluid_params)
            lb_fluid_2 = espressomd.lb.LBFluidGPU(**lb_fluid_params)
        else:
            lb_fluid_1 = espressomd.lb.LBFluid(**lb_fluid_params)
            lb_fluid_2 = espressomd.lb.LBFluid(**lb_fluid_params)

        system.actors.add(lb_fluid_1)
        system.thermostat.set_lb(LB_fluid=lb_fluid_1, gamma=friction_1)

        system.integrator.run(1)

        force_on_part = -friction_1 * np.copy(p.v)

        np.testing.assert_allclose(np.copy(p.f), force_on_part)

        system.integrator.run(100)
        self.assertNotAlmostEqual(lb_fluid_1[3, 3, 3].velocity[0], 0.0)

        system.actors.remove(lb_fluid_1)

        p.v = [1, 0, 0]
        system.integrator.run(0)

        np.testing.assert_allclose(np.copy(p.f), 0.0)

        system.actors.add(lb_fluid_2)
        system.thermostat.set_lb(LB_fluid=lb_fluid_2, gamma=friction_2)

        for pid in itertools.product(range(5), repeat=3):
            np.testing.assert_allclose(
                np.copy(lb_fluid_2[pid].velocity), np.zeros((3,)))

        p.v = [1, 0, 0]

        system.integrator.run(1)

        np.testing.assert_allclose(
            np.copy(p.f), [-friction_2, 0.0, 0.0])

    def test_CPU_LB(self):
        self.switch_test()

    @utx.skipIfMissingGPU()
    def test_GPU_LB(self):
        self.switch_test(GPU=True)


if __name__ == "__main__":
    ut.main()
