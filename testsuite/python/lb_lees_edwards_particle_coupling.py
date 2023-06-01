#
# Copyright (C) 2013-2023 The ESPResSo project
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

system = espressomd.System(box_l=[10, 10, 10])


@utx.skipIfMissingFeatures("WALBERLA")
class LBLeesEdwardsParticleCoupling(ut.TestCase):
    def test_viscous_coupling_with_offset(self):
        system.actors.clear()
        system.time_step = 1
        system.cell_system.skin = 0.1
        system.cell_system.set_n_square()

        offset = 0.9  # (np.random.random()-1/2) * 5*system.box_l[0]
        protocol = lees_edwards.LinearShear(
            shear_velocity=0., initial_pos_offset=offset, time_0=0.)
        system.lees_edwards.set_boundary_conditions(
            shear_direction="x", shear_plane_normal="y", protocol=protocol)

        lbf = espressomd.lb.LBFluidWalberla(
            agrid=1., density=1., kinematic_viscosity=1., tau=system.time_step)
        system.actors.add(lbf)
        system.thermostat.set_lb(LB_fluid=lbf, seed=123, gamma=1)

        idx = int(offset % system.box_l[0])  # lb x index incl offset 
        # the particle is placed in the center of 8 lb points, i.e., 
        # 0.5 away from the planes with the grid points
        pos = [system.box_l[0] / 2., 0., system.box_l[0] / 2.]
        p = system.part.add(pos=pos)
        mid_x = lbf.shape[0] // 2  # lb index for the origin particle
        mid_z = lbf.shape[2] // 2

        upper_y = lbf.shape[1] - 1  # y is shear plane noremal
        # LB Nodes surrounding the particle.
        # On the particle's side 
        unshifted = [lbf[mid_x - 1, 0, mid_z - 1], 
                     lbf[mid_x, 0, mid_z - 1],
                     lbf[mid_x - 1, 0, mid_z],
                     lbf[mid_x, 0, mid_z]]
        # across the Lees Edwads boundary conditoin
        # if offset ^1<.5, it couples to the left neighbor
        # otherwise to the right
        if (offset % 1) >= .5: extra = 1
        else: extra = 0
        shifted_left = [lbf[(mid_x - 1 + idx + extra) % lbf.shape[0], upper_y, mid_z - 1],
                        lbf[(mid_x - 1 + idx + extra) % lbf.shape[0], upper_y, mid_z]]
        shifted_right = [lbf[(mid_x + idx + extra) % lbf.shape[0], upper_y, mid_z - 1],
                         lbf[(mid_x + idx + extra) % lbf.shape[0], upper_y, mid_z]]
        nodes = shifted_left + shifted_right + unshifted

        v0 = np.array([1, 2, 3])
        for n in nodes:
            n.velocity = v0

        system.integrator.run(1)

        # Gather forces applied to the LB by the particle coupling
        lb_forces = np.array([n.last_applied_force for n in nodes])
        lb_force = np.sum(lb_forces, axis=0)

        # total force on lb = - force on particle?
        np.testing.assert_allclose(lb_force, -np.copy(p.f))
        # Particle couples to 8 nodes. On the sie of the particle
        # each lb node should have 1/8 of the force
        for n in unshifted:
            np.testing.assert_allclose(
                np.copy(n.last_applied_force), -np.copy(p.f) / 8)

        # Across the lees edwards boundary, forces ahead and behind
        # the particle in shear directiion are differently distributed
        # For offset 0 (in the middle between two nodes)
        # left and right weighs are equal (0.5)
        weight_r = (offset + .5 - idx) % 1
        weight_l = 1 - weight_r
        np.testing.assert_allclose(weight_l + weight_r, 1)
        for w in (weight_l, weight_r):
            assert w >= 0 and w <= 1
        # The total weight is the product of the weights
        # in all 3 Cartesian directions. These are
        # weight_l/r in the shear direction and 1/2 in the other two
        for n in shifted_left:
            np.testing.assert_allclose(
                np.copy(n.last_applied_force), -np.copy(p.f) * 1 / 2 * 1 / 2 * weight_l)
        for n in shifted_right:
            np.testing.assert_allclose(
                np.copy(n.last_applied_force), -np.copy(p.f) * 1 / 2 * 1 / 2 * weight_r)

        # Check the LB velocity interpolation
        # at a positoin, where the ghost layer of the
        # lees edwards shear plane contributes
        lbf[:, :, :].velocity = np.zeros(3)
        lbf[:, :, :].velocity = [0, 0, 0]

        lower_nodes = unshifted
        upper_nodes = shifted_left + shifted_right
        for n in lower_nodes:
            n.velocity = v0
        for n in upper_nodes:
            n.velocity = - v0
        # When the offset modulo 1 is not zero
        # in the ghost layer, neighboring cell swith zero velocity contribute
        expected_v = 1 / 2 * v0 +\
            1 / 4 * -v0 + 1 / 4 * (1 - abs(1 - (offset % 1))) * -v0
        p.update(dict(pos=pos, v=np.zeros(3)))
        np.testing.assert_allclose(
            np.copy(lbf.get_interpolated_velocity(pos=pos)),
            expected_v)

    def test_viscous_coupling_with_shear_vel(self):
        # Places a co-moving particle close to the LE boundary
        # in shear flow. chesk that it remains force free
        # this is only the case, if the periodic imagesin the 
        # halo regoin calculate the drag force including the LE
        # shear velocity.
        system.actors.clear()
        system.time_step = 0.1
        system.cell_system.skin = 0.1
        system.cell_system.set_n_square()
        v_shear = 0.5
        protocol = lees_edwards.LinearShear(
            shear_velocity=v_shear, initial_pos_offset=0, time_0=0.)
        system.lees_edwards.set_boundary_conditions(
            shear_direction="x", shear_plane_normal="y", protocol=protocol)

        lbf = espressomd.lb.LBFluidWalberla(
            agrid=1., density=1., kinematic_viscosity=1., tau=system.time_step)
        system.actors.add(lbf)
        system.thermostat.set_lb(LB_fluid=lbf, seed=123, gamma=1)
        system.integrator.run(200) 
        p = system.part.add(
            pos=(
                0, 0, 0), v=lbf.get_interpolated_velocity(
                pos=(
                    0, 0, 0)))
        for _ in range(100): 
            system.integrator.run(1)
            np.testing.assert_allclose(p.f, np.zeros(3), atol=1E-7)


if __name__ == '__main__':
    ut.main()
