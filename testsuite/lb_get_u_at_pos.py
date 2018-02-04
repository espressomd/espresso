import sys
import unittest as ut
import numpy as np
import numpy.testing
import espressomd
from espressomd import lb


@ut.skipIf(not espressomd.has_features("LB_GPU") or espressomd.has_features(
    "SHANCHEN"), "LB_GPU feature not available, skipping test!")
class TestLBGetUAtPos(ut.TestCase):
    """
    Check velocities at particle positions are sorted by ``id`` and
    quantitatively correct (only LB GPU).

    """
    @classmethod
    def setUpClass(self):
        self.params = {
            'tau': 0.01,
            'agrid': 0.5,
            'box_l': [12.0, 12.0, 12.0],
            'dens': 0.85,
            'viscosity': 30.0,
            'friction': 2.0,
            'gamma': 1.5
        }
        self.system = espressomd.System(box_l=[1.0, 1.0, 1.0])
        self.system.seed = self.system.cell_system.get_state()['n_nodes'] * [1234]
        self.system.box_l = self.params['box_l']
        self.system.cell_system.skin = 0.4
        self.system.time_step = 0.01
        self.n_nodes_per_dim = int(self.system.box_l[0] / self.params['agrid'])
        for p in range(self.n_nodes_per_dim):
            # Set particles exactly between two LB nodes in x direction.
            self.system.part.add(id=p,
                                 pos=[(p + 1) * self.params['agrid'],
                                      0.5 * self.params['agrid'],
                                      0.5 * self.params['agrid']])
        self.lb_fluid = lb.LBFluidGPU(
            visc=self.params['viscosity'],
            dens=self.params['dens'],
            agrid=self.params['agrid'],
            tau=self.params['tau'],
            fric=self.params['friction']
        )
        self.system.actors.add(self.lb_fluid)
        self.vels = np.zeros((self.n_nodes_per_dim, 3))
        self.vels[:, 0] = np.arange(self.n_nodes_per_dim, dtype=float)
        self.interpolated_vels = self.vels.copy()
        self.interpolated_vels[:, 0] += 0.5
        for n in range(self.n_nodes_per_dim):
            self.lb_fluid[n, 0, 0].velocity = self.vels[n, :]
        self.system.integrator.run(0)

    def test_get_u_at_pos(self):
        """
        Test if linear interpolated velocities are equal to the velocities at
        the particle positions. This test uses the two-point coupling under
        the hood.

        """
        numpy.testing.assert_allclose(
                self.interpolated_vels[:-1],
            self.lb_fluid.get_interpolated_fluid_velocity_at_positions(
                self.system.part[:].pos)[:-1],
            atol=1e-4)


if __name__ == "__main__":
    suite = ut.TestSuite()
    suite.addTests(ut.TestLoader().loadTestsFromTestCase(TestLBGetUAtPos))
    result = ut.TextTestRunner(verbosity=4).run(suite)
    sys.exit(not result.wasSuccessful())
