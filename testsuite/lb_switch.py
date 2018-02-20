from __future__ import print_function
import unittest as ut
import numpy as np
import espressomd
import espressomd.lb
from espressomd import *
from tests_common import abspath
from itertools import product

@ut.skipIf(not espressomd.has_features(["LB"]),
           "Features not available, skipping test!")
class LBSwitchActor(ut.TestCase):
    system = espressomd.System(box_l=[10.0, 10.0, 10.0])
    n_nodes = system.cell_system.get_state()["n_nodes"]
    system.seed = range(n_nodes)

    system.time_step = 0.01
    system.cell_system.skin = 0.1

    def test(self):
        system = self.system

        system.part.add(pos=[5.,5.,5.], v=[1.,0,0])
        gamma_1 = 1.0
        gamma_2 = 2.0

        system.thermostat.set_lb(kT=0.0)

        lb_fluid_1 = espressomd.lb.LBFluid(
            agrid=2.0, dens=1.0, visc=1.0, fric=gamma_1, tau=0.03)

        lb_fluid_2 = espressomd.lb.LBFluid(
            agrid=2.0, dens=1.0, visc=1.0, fric=gamma_2, tau=0.03)

        system.actors.add(lb_fluid_1)

        system.integrator.run(1)

        np.testing.assert_allclose(system.part[0].f, [-gamma_1, 0.0, 0.0])

        system.integrator.run(100)
        self.assertNotAlmostEqual(lb_fluid_1[3,3,3].velocity[0], 0.0)

        system.actors.remove(lb_fluid_1)

        system.part[0].v = [1,0,0]
        system.integrator.run(0)

        np.testing.assert_allclose(system.part[0].f, 0.0)

        system.actors.add(lb_fluid_2)
        for p in product(range(5), range(5), range(5)):
            np.testing.assert_allclose(lb_fluid_2[p].velocity, np.zeros((3,)))

        system.part[0].v = [1,0,0]

        system.integrator.run(1)

        np.testing.assert_allclose(system.part[0].f, [-gamma_2, 0.0, 0.0])

if __name__ == "__main__":
    ut.main()

