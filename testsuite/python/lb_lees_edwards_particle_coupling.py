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

def min_image_dist(a,b,l):
   res = b-a 
   for i in range(3):
     while res[i]<-l[i]/2: res+=l[i]
     while res[i]>=l[i]/2: res-=l[i]
   return res

def coupling_weight(pos_lb_units, node_idx, lb_shape):
    # 1. For the coupling weights it does not matter on which side of the lb_node the posiiton is
    # 2. To determine the lb node to position distance, we need 
    # minimum image convetion, node and coupling position can be at different sides
    # of a periodic boudnary  
    dx = np.abs(min_image_dist(pos_lb_units,node_idx, lb_shape))
    # If the coupling point is >=1 lattice constant away from the node, no coupling.
    if np.any(dx>=1): 
       weight=0
    else:
      # distance pos to node via module with lattice constant 1
      weight =np.product(1-dx)

    return weight


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

    def check_velocity_interpolation(self, pos_offset, shear_vel, test_positions):
        system.lb = None
        system.time_step = 1
        system.cell_system.skin = 0.1
        system.cell_system.set_n_square()
        protocol = lees_edwards.LinearShear(
            shear_velocity=shear_vel, initial_pos_offset=pos_offset, time_0=0.)
        system.lees_edwards.set_boundary_conditions(
            shear_direction="x", shear_plane_normal="y", protocol=protocol)
        lbf = espressomd.lb.LBFluidWalberla(
            agrid=1., density=1., kinematic_viscosity=1., tau=system.time_step)
        system.lb = lbf
        system.thermostat.set_lb(LB_fluid=lbf, seed=123, gamma=1)
        system.part.clear()
        v_x = lambda x: np.interp(x,[0.5,lbf.shape[0]-.5],[0,lbf.shape[0]-1],period=lbf.shape[0])
        nodes_at_y_boundary = list(lbf[:,0,:]) +list(lbf[:,lbf.shape[1]-1,:])
        for n in nodes_at_y_boundary:
                node_x = 0.5+n.index[0]
                n.velocity = [v_x(node_x), 0, 0]
        for pos in test_positions:

            y = pos[1]
            if abs(y <= 0.5): 
               pref = -1
               dist_to_unshifted_lb_nodes = 0.5-y
            elif y >= system.box_l[1] - .5: 
               pref = 1
               dist_to_unshifted_lb_nodes = y-(system.box_l[2]-.5)
            else: raise Exception()
            vel_shift = pref * shear_vel
            xs = 0.5+np.arange(lbf.shape[0])
            ys = [v_x(x-pref*pos_offset) for x in xs]
            v_x_shifted = lambda x: np.interp(x,xs,ys,period=system.box_l[0])
            unshifted_vel = v_x(pos[0])
            shifted_vel = v_x_shifted(pos[0]) +vel_shift
            weight_unshifted = 1-dist_to_unshifted_lb_nodes
            weight_shifted = 1-weight_unshifted
            expected_vel = np.array([weight_unshifted *unshifted_vel +weight_shifted * shifted_vel,0,0])
            observed_vel = np.copy(lbf.get_interpolated_velocity(pos=pos))
            np.testing.assert_allclose(observed_vel, expected_vel)
    
    def test_vel_interpol_all(self):
       n = 25
       xs = np.linspace(0,system.box_l[0],n)
       y_ls = [0.2]*n
       y_us = [system.box_l[1]-.2]*n
       zs = np.random.random(n) * system.box_l[2]
       pos_lower = np.vstack((xs,y_ls,zs)).T       
       pos_upper = np.vstack((xs,y_us,zs)).T       
       pos_all = np.vstack((pos_lower, pos_upper))
       # non-integer offset 
       pos_offsets = 100 * system.box_l[0] *(np.random.random(10) -.5)
       for pos_offset in pos_offsets:
         self.check_velocity_interpolation(pos_offset, 0, pos_all)

    def test_viscous_coupling_with_shear_vel(self):
        # Places a co-moving particle close to the LE boundary
        # in shear flow. checks that it remains force free
        # this is only the case, if the periodic images in the 
        # halo regoin calculate the drag force including the LE
        # shear velocity.
        system.lb = None
        system.part.clear()
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
        system.integrator.run(5000) 
        for n in lbf[:,:,:]:
           np.testing.assert_allclose(n.velocity[1:],[0,0],atol=1E-8)
        pos = np.random.random(3) * system.box_l
        p = system.part.add(
            pos=pos, v=lbf.get_interpolated_velocity(pos=pos))
        np.testing.assert_allclose(p.v[1:],[0,0],atol=1E-8)
        for _ in range(1000): 
            system.integrator.run(1)
            np.testing.assert_allclose(np.copy(p.f), np.zeros(3), atol=2E-6)

    def test_momentum_conservation(self):
        system.lb = None
        system.part.clear()
        system.time_step = 0.01
        system.cell_system.skin = 0.1
        system.cell_system.set_n_square()
        v_shear = 0.5
        protocol = lees_edwards.LinearShear(
            shear_velocity=0.5, initial_pos_offset=13.7, time_0=0.)
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
        for i in range(100): 
            before = (p.pos_folded,p.v,lbf.get_interpolated_velocity(pos=p.pos_folded))
            system.integrator.run(1)
            after = (p.pos_folded,p.v,lbf.get_interpolated_velocity(pos=p.pos_folded))
            np.testing.assert_allclose(-np.copy(p.f), np.copy(np.sum(lbf[:,:,:].last_applied_force,axis=(0,1,2))),atol=1E-9)
            print("b", before)
            print("a",after)
            current_mom = np.copy(system.analysis.linear_momentum())
            print("m" ,(initial_mom-current_mom)[1:])
            print()
            np.testing.assert_allclose(initial_mom[1:], current_mom[1:], atol=2E-7)


if __name__ == '__main__':
    ut.main()
