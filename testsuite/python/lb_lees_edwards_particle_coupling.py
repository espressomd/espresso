#
# Copyright (C) 2013-2019 The ESPResSo project
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
import espressomd.lees_edwards as lees_edwards
import espressomd
import espressomd.lb
import numpy as np
import unittest_decorators as utx


@utx.skipIfMissingFeatures("WALBERLA")
class LBLeesEdwardsParticleCoupling(ut.TestCase):
    def test(self):
        system = espressomd.System(box_l=[10, 10, 10])

        system.time_step = 1
        system.cell_system.skin = 0.1
        system.cell_system.set_n_square()

        offset = 1
        idx = int(offset)
        protocol = lees_edwards.LinearShear(
            shear_velocity=0., initial_pos_offset=offset, time_0=0.)
        system.lees_edwards.set_boundary_conditions(
            shear_direction="x", shear_plane_normal="y", protocol=protocol)

        lbf = espressomd.lb.LBFluidWalberla(
            agrid=1, density=1, viscosity=1, tau=system.time_step)
        system.actors.add(lbf)
        system.thermostat.set_lb(LB_fluid=lbf, seed=123, gamma=1)

        pos = [system.box_l[0] / 2., 0., system.box_l[0] / 2.]
        p = system.part.add(pos=pos)
        v0 = np.array([1, 2, 3])
        mid_x = lbf.shape[0] // 2
        mid_z = lbf.shape[2] // 2

        upper_y = lbf.shape[1] - 1
        nodes = [lbf[mid_x - 1, 0, mid_z],
                 lbf[mid_x, 0, mid_z - 1],
                 lbf[mid_x - 1, 0, mid_z],
                 lbf[mid_x, 0, mid_z],
                 lbf[mid_x - 1 + idx, upper_y, mid_z],
                 lbf[mid_x + idx, upper_y, mid_z - 1],
                 lbf[mid_x - 1 + idx, upper_y, mid_z],
                 lbf[mid_x + idx, upper_y, mid_z]]
        for n in nodes:
            n.velocity = v0

        system.integrator.run(1)
        lb_forces = np.array([n.last_applied_force for n in nodes])
        lb_force = np.sum(lb_forces, axis=0)
        np.testing.assert_allclose(lb_force, -np.copy(p.f))
        for f in lb_forces:
            np.testing.assert_allclose(f, lb_forces[0])

        lbf[:, :, :].velocity = [0, 0, 0]

        lower_nodes = nodes[:4]
        upper_nodes = nodes[4:]
        for n in lower_nodes:
            n.velocity = v0
        for n in upper_nodes:
            n.velocity = - v0
        p.update(dict(pos=pos, v=np.zeros(3)))
        np.testing.assert_allclose(
            np.copy(lbf.get_interpolated_velocity(pos)),
            np.zeros(3))
        system.integrator.run(1)
        np.testing.assert_allclose(np.copy(p.pos), pos)
        np.testing.assert_allclose(np.copy(p.f), np.zeros(3))
        for n in nodes:
            np.testing.assert_allclose(
                np.copy(n.last_applied_force), np.zeros(3))


if __name__ == '__main__':
    ut.main()
