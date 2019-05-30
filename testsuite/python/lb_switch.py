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
from __future__ import print_function
import unittest as ut
import numpy as np
import espressomd
import espressomd.lb
from tests_common import abspath
from itertools import product


@ut.skipIf(not espressomd.has_features(["EXTERNAL_FORCES"]),
           "Features not available, skipping test!")
class LBSwitchActor(ut.TestCase):
    system = espressomd.System(box_l=[10.0, 10.0, 10.0])

    system.time_step = 0.01
    system.cell_system.skin = 0.1

    def switch_test(self, GPU=False):
        system = self.system
        system.actors.clear()
        system.part.add(pos=[1., 1., 1.], v=[1., 0, 0], fix=[1, 1, 1])
        ext_force_density = [0.2, 0.3, 0.15]

        lb_fluid_params = {
            'agrid': 2.0, 'dens': 1.0, 'visc': 1.0, 'tau': 0.03}
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

        force_on_part = -friction_1 * np.copy(system.part[0].v)

        np.testing.assert_allclose(np.copy(system.part[0].f), force_on_part)

        system.integrator.run(100)
        self.assertNotAlmostEqual(lb_fluid_1[3, 3, 3].velocity[0], 0.0)

        system.actors.remove(lb_fluid_1)

        system.part[0].v = [1, 0, 0]
        system.integrator.run(0)

        np.testing.assert_allclose(np.copy(system.part[0].f), 0.0)

        system.actors.add(lb_fluid_2)
        system.thermostat.set_lb(LB_fluid=lb_fluid_2, gamma=friction_2)

        for p in product(range(5), range(5), range(5)):
            np.testing.assert_allclose(
                np.copy(lb_fluid_2[p].velocity), np.zeros((3,)))

        system.part[0].v = [1, 0, 0]

        system.integrator.run(1)

        np.testing.assert_allclose(
            np.copy(system.part[0].f), [-friction_2, 0.0, 0.0])

    def test_CPU_LB(self):
        self.switch_test()

    @ut.skipIf((not espressomd.gpu_available() or not espressomd.has_features(["CUDA"])
                ),
               "CUDA not available or no gpu present, skipping test.")
    def test_GPU_LB(self):
        self.switch_test(GPU=True)


if __name__ == "__main__":
    ut.main()
