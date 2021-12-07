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
import sys
import unittest as ut
import unittest_decorators as utx
import numpy as np
import numpy.testing
import espressomd
import espressomd.lb


@utx.skipIfMissingGPU()
class TestLBGetUAtPos(ut.TestCase):

    """
    Check velocities at particle positions are sorted by ``id`` and
    quantitatively correct (only LB GPU).

    """
    @classmethod
    def setUpClass(cls):
        cls.params = {
            'tau': 0.01,
            'agrid': 0.5,
            'box_l': [12.0, 12.0, 12.0],
            'dens': 0.85,
            'viscosity': 30.0,
            'friction': 2.0,
            'gamma': 1.5
        }
        cls.system = espressomd.System(box_l=[1.0, 1.0, 1.0])
        cls.system.box_l = cls.params['box_l']
        cls.system.cell_system.skin = 0.4
        cls.system.time_step = 0.01
        cls.n_nodes_per_dim = int(cls.system.box_l[0] / cls.params['agrid'])
        for p in range(cls.n_nodes_per_dim):
            # Set particles exactly between two LB nodes in x direction.
            cls.system.part.add(pos=[(p + 1) * cls.params['agrid'],
                                     0.5 * cls.params['agrid'],
                                     0.5 * cls.params['agrid']])
        cls.lb_fluid = espressomd.lb.LBFluidGPU(
            visc=cls.params['viscosity'],
            dens=cls.params['dens'],
            agrid=cls.params['agrid'],
            tau=cls.params['tau'],
        )
        cls.system.actors.add(cls.lb_fluid)
        cls.vels = np.zeros((cls.n_nodes_per_dim, 3))
        cls.vels[:, 0] = np.arange(cls.n_nodes_per_dim, dtype=float)
        cls.interpolated_vels = cls.vels.copy()
        cls.interpolated_vels[:, 0] += 0.5
        for n in range(cls.n_nodes_per_dim):
            cls.lb_fluid[n, 0, 0].velocity = cls.vels[n, :]
        cls.system.integrator.run(0)

    def test_get_u_at_pos(self):
        """
        Test if linear interpolated velocities are equal to the velocities at
        the particle positions. This test uses the two-point coupling under
        the hood.

        """
        numpy.testing.assert_allclose(
            self.interpolated_vels[:-1],
            self.lb_fluid.get_interpolated_fluid_velocity_at_positions(
                self.system.part.all().pos, False)[:-1],
            atol=1e-4)


if __name__ == "__main__":
    suite = ut.TestSuite()
    suite.addTests(ut.TestLoader().loadTestsFromTestCase(TestLBGetUAtPos))
    result = ut.TextTestRunner(verbosity=4).run(suite)
    sys.exit(not result.wasSuccessful())
