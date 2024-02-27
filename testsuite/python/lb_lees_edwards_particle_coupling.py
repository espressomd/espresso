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
import itertools
import copy
import unittest_decorators as utx
from tests_common import fold_index

system = espressomd.System(box_l=[10, 10, 10])


def unit_vec(k):
    res = np.zeros(3)
    res[k] = 1
    return res


def within_grid(idx, shape):
    return np.all(idx >= 0) and np.all(idx < shape)


def coupling_weight(pos_lb_units, node_idx, lb_shape):
    # minimum image convention
    dx = np.abs(pos_lb_units - node_idx) 
    for i in range(3):
        while dx[i] > lb_shape[i] / 2: dx[i] -= lb_shape[i]

    if np.any(np.abs(dx) >= 1): return 0.  # only couple to neighbors

    return np.product(1 - np.abs(dx))


class MockLBF:
    shape = np.array((10, 10, 10))
    agrid = 1

    def __getitem__(self, idx):
        assert within_grid(idx, self.shape)
        return f"node {idx}"


mock_lbf = MockLBF()


def le_aware_lb_nodes_around_pos(
        folded_pos, lbf, le_pos_offset, shear_direction, shear_plane_normal):
    """Returns LB node(s) relevant for interpolation around the given position"""

    # helper to get lb node from index with folding
    def lb_node(idx): return lbf[fold_index(idx, lbf.shape)]

    # center of first lb lattice point is at 0.5 in each Cartesian diereciotn.

    # determine position in lb units for 3 cases;
    # unshifted, lees-edwards shifted to the left and to the right
    pos_unshifted_lb_units = folded_pos / lbf.agrid - 0.5  # relative to node centers
    shear_vec = unit_vec(shear_direction)
    pos_shifted_left_lb_units = (
        folded_pos - shear_vec * le_pos_offset) / lbf.agrid - 0.5
    pos_shifted_right_lb_units = (
        folded_pos + shear_vec * le_pos_offset) / lbf.agrid - 0.5

    # Particle couples to its 8 neighboring lattice sites
    # when a particle is at a lattice site in any coordinate, the right hand side neighbor is included
    # but the coupling weight for that neighbor is 0

    # find the lower left lb node to which the particle couples
    lower_idx_unshifted = np.array(np.floor(pos_unshifted_lb_units), dtype=int)
    lower_idx_shifted_left = np.array(
        np.floor(pos_shifted_left_lb_units), dtype=int)
    lower_idx_shifted_right = np.array(
        np.floor(pos_shifted_right_lb_units), dtype=int)

    ijks = np.array(list(itertools.product([0, 1], repeat=3)))
    indices_unshifted = [lower_idx_unshifted + ijk for ijk in ijks]
    indices_shifted_left = [lower_idx_shifted_left + ijk for ijk in ijks]
    indices_shifted_right = [lower_idx_shifted_right + ijk for ijk in ijks]

    # Nodes with an index within the primary box in shear_plane_normal direction
    # do not need Lees-Edwards handling
    dont_need_shift = [
        idx for idx in indices_unshifted if within_grid(
            idx[shear_plane_normal],
            lbf.shape[shear_plane_normal])]
    unshifted_nodes = [lb_node(idx) for idx in dont_need_shift]
    unshifted_weights = [
        coupling_weight(
            pos_unshifted_lb_units,
            idx,
            lbf.shape) for idx in dont_need_shift]

    # Handle indices which are not in the primary box in the sheare plane
    # normal
    to_be_shifted_left = [
        (idx,
         ijk) for idx,
        ijk in zip(
            indices_unshifted,
            ijks) if idx[shear_plane_normal] >= lbf.shape[shear_plane_normal]]
    to_be_shifted_right = [(idx, ijk) for idx, ijk in zip(
        indices_unshifted, ijks) if idx[shear_plane_normal] < 0]

    # replace the index in shear direction   
    shifted_left = copy.deepcopy(to_be_shifted_left)    
    for idx, ijk in shifted_left:
        idx[shear_direction] = lower_idx_shifted_left[shear_direction] + \
            ijk[shear_direction]
    shifted_right = copy.deepcopy(to_be_shifted_right)    
    for idx, ijk in shifted_right:
        idx[shear_direction] = lower_idx_shifted_right[shear_direction] + \
            ijk[shear_direction]

    weights_shifted_left = [
        coupling_weight(
            pos_shifted_left_lb_units,
            idx,
            lbf.shape) for idx,
        _ in shifted_left]
    weights_shifted_right = [
        coupling_weight(
            pos_shifted_right_lb_units,
            idx,
            lbf.shape) for idx,
        _ in shifted_right]

    shifted_nodes = [lb_node(idx) for idx, _ in (shifted_left + shifted_right)]
    shifted_weights = weights_shifted_left + weights_shifted_right

    np.testing.assert_allclose(sum(unshifted_weights + shifted_weights), 1)
    assert len(shifted_nodes + unshifted_nodes) == 8
    assert len(set(shifted_nodes + unshifted_nodes)) == 8  # no duplicates

    return unshifted_nodes, shifted_nodes, unshifted_weights, shifted_weights

# print(le_aware_lb_nodes_around_pos(np.array((0.5,9.99999,0)),mock_lbf,0.51,0,1))
# exit()


@utx.skipIfMissingFeatures("WALBERLA")
class LBLeesEdwardsParticleCoupling(ut.TestCase):
    def test_viscous_coupling_with_offset(self):
        system.lb = None
        system.time_step = 1
        system.cell_system.skin = 0.1
        system.cell_system.set_n_square()
        offset = 2 * (np.random.random() - 1) * 3 * system.box_l[1]
        protocol = lees_edwards.LinearShear(
            shear_velocity=0., initial_pos_offset=offset, time_0=0.)
        system.lees_edwards.set_boundary_conditions(
            shear_direction="x", shear_plane_normal="y", protocol=protocol)
        lbf = espressomd.lb.LBFluidWalberla(
            agrid=1., density=1., kinematic_viscosity=1., tau=system.time_step)
        system.lb = lbf
        system.thermostat.set_lb(LB_fluid=lbf, seed=123, gamma=1)
        for _ in range(10):
            system.part.clear()
            lbf[:, :, :].velocity = np.zeros(3)

            x = np.random.random() * system.box_l[0]
            z = np.random.random() * system.box_l[2]
            # within 0.5 of the le boundar       pos = np.array(x,y,z)
            y = (np.random.random() - .5) % system.box_l[1]
            pos = np.array((x, y, z))
            p = system.part.add(pos=pos)
            v0 = np.random.random(3) - 1 / 2

            nodes_unshifted, nodes_shifted, weights_unshifted, weights_shifted = \
                le_aware_lb_nodes_around_pos(pos, lbf, offset, 0, 1)
            all_nodes = nodes_unshifted + nodes_shifted
            all_weights = weights_unshifted + weights_shifted

            for n in all_nodes:
                n.velocity = v0

            system.integrator.run(1)

            # validate our assumptions about which lb nodes get a force
            # from the coupling. Eactly the nodes listed in ` nodes` 
            # should have received a force during coupling.
            lb_nodes_with_force_idx = sorted(
                [n.index for n in lbf[:, :, :] if np.any(n.last_applied_force != 0)])
            expected_nodes_idx = sorted(
                [n.index for n, w in zip(all_nodes, all_weights) if w > 0])
            np.testing.assert_array_equal(
                lb_nodes_with_force_idx, expected_nodes_idx)

            # Gather forces applied to the LB by the particle coupling
            lb_force = np.sum(
                np.array([n.last_applied_force for n in all_nodes]), axis=0)

            # total force on lb = - force on particle?
            np.testing.assert_allclose(lb_force, -np.copy(p.f))

            # force on individual nodes
            for n, w in zip(all_nodes, all_weights):
                np.testing.assert_allclose(
                    np.copy(n.last_applied_force), -w * np.copy(p.f))

    def test_velocity_interpolation(self):
        system.lb = None
        system.time_step = 1
        system.cell_system.skin = 0.1
        system.cell_system.set_n_square()
        # 2 * (np.random.random() - 1) * 3 * system.box_l[1]
        pos_offset = 1.6 
        shear_vel = 0  # np.random.random()-1/2
        protocol = lees_edwards.LinearShear(
            shear_velocity=shear_vel, initial_pos_offset=pos_offset, time_0=0.)
        system.lees_edwards.set_boundary_conditions(
            shear_direction="x", shear_plane_normal="y", protocol=protocol)
        lbf = espressomd.lb.LBFluidWalberla(
            agrid=1., density=1., kinematic_viscosity=1., tau=system.time_step)
        system.lb = lbf
        system.thermostat.set_lb(LB_fluid=lbf, seed=123, gamma=1)
        system.part.clear()
        for x in np.linspace(0, system.box_l[0], 21):
            for n in lbf[:, :, :]:
                n.velocity = [n.index[0], 0, 0]

            z = 4  # np.random.random()*system.box_l[2]
            y = 0
            pos = np.array((x, y, z))
            nodes_unshifted, nodes_shifted, weights_unshifted, weights_shifted = \
                le_aware_lb_nodes_around_pos(pos, lbf, pos_offset, 0, 1)
            all_nodes = nodes_unshifted + nodes_shifted
            # all_weights = weights_unshifted + weights_shifted

            vels_unshifted_nodes = np.array(
                [n.velocity for n in nodes_unshifted])
            vels_shifted_nodes = np.array([n.velocity for n in nodes_shifted])
#            for n, v, w in zip(
#                    nodes_unshifted, vels_unshifted_nodes, weights_unshifted):
#                print(n.index, v, w, v * w)
#            for n, v, w in zip(
#                    nodes_shifted, vels_shifted_nodes, weights_shifted):
#                print(n.index, v, w, v * w)

            if abs(y <= 0.5): vels_shifted_nodes[:, 0] -= shear_vel 
            elif y >= system.box_l[1] - .5: vels_shifted_nodes[:, 0] += shear_vel 
            # else: raise Exception()

            vel_contrib_unshifted = np.array(
                [w * vel for vel, w in zip(vels_unshifted_nodes, weights_unshifted)])
            vel_contrib_shifted = np.array(
                [w * vel for vel, w in zip(vels_shifted_nodes, weights_shifted)])
            expected_vel = np.sum(vel_contrib_unshifted,
                                  axis=0) + np.sum(vel_contrib_shifted, axis=0)

            observed_vel = np.copy(lbf.get_interpolated_velocity(pos=pos))
            print(pos[0], pos_offset, observed_vel, expected_vel)
            np.testing.assert_allclose(observed_vel, expected_vel)

    def atest_viscous_coupling_with_shear_vel(self):
        # Places a co-moving particle close to the LE boundary
        # in shear flow. checks that it remains force free
        # this is only the case, if the periodic images in the 
        # halo regoin calculate the drag force including the LE
        # shear velocity.
        system.lb = None
        system.time_step = 0.1
        system.cell_system.skin = 0.1
        system.cell_system.set_n_square()
        v_shear = 2 * (np.random.random() - 1 / 2) 
        protocol = lees_edwards.LinearShear(
            shear_velocity=v_shear, initial_pos_offset=(np.random.random() - 1 / 2) * 5 * system.box_l[0], time_0=np.random.random())
        system.lees_edwards.set_boundary_conditions(
            shear_direction="x", shear_plane_normal="y", protocol=protocol)

        lbf = espressomd.lb.LBFluidWalberla(
            agrid=1., density=1., kinematic_viscosity=1., tau=system.time_step)
        system.lb = lbf
        system.thermostat.set_lb(LB_fluid=lbf, seed=123, gamma=1)
        system.integrator.run(2000) 
        pos = np.random.random(3) * system.box_l
        p = system.part.add(
            pos=pos, v=lbf.get_interpolated_velocity(pos=pos))
        for _ in range(1000): 
            system.integrator.run(1)
            np.testing.assert_allclose(np.copy(p.f), np.zeros(3), atol=2E-6)

    def xtest_momentum_conservation(self):
        system.lb = None
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
        system.lb = lbf
        system.thermostat.set_lb(LB_fluid=lbf, seed=123, gamma=1)
        pos = (0, 0, 0)
        p = system.part.add(
            pos=pos, v=(0, 0, 0))
        system.integrator.run(1) 
        initial_mom = np.copy(system.analysis.linear_momentum())
        for _ in range(1000): 
            system.integrator.run(1)
            current_mom = np.copy(system.analysis.linear_momentum())
            print(current_mom, p.pos_folded)
            np.testing.assert_allclose(initial_mom, current_mom, atol=1E-6)


if __name__ == '__main__':
    ut.main()
