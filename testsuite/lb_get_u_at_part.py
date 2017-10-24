import sys
import unittest as ut
import numpy as np
import numpy.testing
import espressomd
from espressomd import lb


@ut.skipIf(not espressomd.has_features("LB_GPU") or espressomd.has_features("SHANCHEN"),
           "LB_GPU feature not available, skipping test!")
class TestLBGetUAtPart(ut.TestCase):
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
        self.system = espressomd.System()
        self.system.box_l = self.params['box_l']
        self.system.cell_system.skin = 0.4
        self.system.time_step = 0.01
        self.n_nodes_per_dim = int(self.system.box_l[0] / self.params['agrid'])
        for p in range(self.n_nodes_per_dim):
            if p % 2 == 0:
                self.system.part.add(id=p,
                                     pos=[0.5 * self.params['agrid'] * (p + 1),
                                          0.5 * self.params['agrid'],
                                          0.5 * self.params['agrid']])
        self.lb_fluid = lb.LBFluid_GPU(
            visc=self.params['viscosity'],
            dens=self.params['dens'],
            agrid=self.params['agrid'],
            tau=self.params['tau'],
            fric=self.params['friction']
        )
        self.system.actors.add(self.lb_fluid)
        self.vels = np.zeros((self.n_nodes_per_dim, 3))
        self.vels[:, 0] = np.arange(self.n_nodes_per_dim, dtype=float)
        for n in range(self.n_nodes_per_dim):
            self.lb_fluid[n, 0, 0].velocity = self.vels[n, :]
        self.system.integrator.run(0)

    def test_get_u_at_part_two_point(self):
        """
        Test if node velocities are equal to the velocities at the particle positions.
        This test uses the two-point coupling under the hood.

        """
        numpy.testing.assert_allclose(
            self.vels[:self.vels.shape[0] / 2, :], self.lb_fluid.get_fluid_velocity_at_particle_positions(
                particle_coupling="2pt"), atol=1e-6)

    #def test_get_u_at_part_three_point(self):
    #    """
    #    Test if node velocities are equal to the velocities at the particle positions.
    #    This test uses the three-point coupling under the hood.

    #    """
    #    numpy.testing.assert_allclose(
    #        self.vels[:self.vels.shape[0] / 2, :], self.lb_fluid.get_fluid_velocity_at_particle_positions(
    #            particle_coupling="3pt"), atol=1e-6)


if __name__ == "__main__":
    suite = ut.TestSuite()
    suite.addTests(ut.TestLoader().loadTestsFromTestCase(TestLBGetUAtPart))
    result = ut.TextTestRunner(verbosity=4).run(suite)
    sys.exit(not result.wasSuccessful())
