import sys
import numpy as np
import unittest as ut

import espressomd
import espressomd.lb


@ut.skipIf(
    not espressomd.has_features(
        ['LB_GPU']) or espressomd.has_features("SHANCHEN"),
           "Features not available, skipping test!")
class LBGPUViscous(ut.TestCase):
    system = espressomd.System(box_l=[10.0] * 3)
    system.time_step = 0.01
    system.cell_system.skin = 0.4
    agrid = 0.5
    dens = 0.85
    viscosity = 30.0
    friction = 3.5

    def test_viscous_coupling(self):
        self.system.thermostat.turn_off()
        self.system.actors.clear()
        self.system.part.clear()
        v_part = np.array([1, 2, 3])
        v_fluid = np.array([1.2, 4.3, 0.2])
        self.lbf = espressomd.lb.LBFluidGPU(
            visc=self.viscosity, dens=self.dens, agrid=self.agrid, tau=self.system.time_step, fric=self.friction)
        self.system.actors.add(self.lbf)
        self.system.part.add(
            pos=[0.5 * self.agrid] * 3, v=v_part, fix=[1, 1, 1])
        self.lbf[0, 0, 0].velocity = v_fluid
        self.system.integrator.run(1)
        np.testing.assert_allclose(
            np.copy(self.system.part[0].f), -self.friction * (v_part - v_fluid), atol=1e-3)


if __name__ == "__main__":
    suite = ut.TestSuite()
    suite.addTests(ut.TestLoader().loadTestsFromTestCase(LBGPUViscous))
    result = ut.TextTestRunner(verbosity=4).run(suite)
    sys.exit(not result.wasSuccessful())
