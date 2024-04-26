#
# Copyright (C) 2013-2022 The ESPResSo project
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
import numpy as np

import espressomd
import espressomd.lb
import espressomd.shapes
import espressomd.utils
import espressomd.propagation

import tests_common
import unittest_decorators as utx


class VirtualSitesTracersCommon:
    agrid = 0.5
    box_height = 10. * agrid
    box_lw = 8. * agrid
    system = espressomd.System(box_l=(box_lw, box_lw, box_height))
    system.time_step = 0.08
    system.cell_system.skin = 0.1

    def setUp(self):
        self.system.box_l = (self.box_lw, self.box_lw, self.box_height)

    def tearDown(self):
        self.system.lb = None
        self.system.part.clear()
        self.system.thermostat.turn_off()

    def set_lb(self, ext_force_density=(0, 0, 0), dir_walls=2):
        self.system.lb = None
        self.lbf = self.LBClass(
            kT=0.0, agrid=self.agrid, density=1., kinematic_viscosity=1.8,
            tau=self.system.time_step, ext_force_density=ext_force_density)
        self.system.lb = self.lbf
        self.system.thermostat.set_lb(LB_fluid=self.lbf, gamma=1.)

        # Setup boundaries
        normal = [0, 0, 0]
        normal[dir_walls] = 1
        wall_shape = espressomd.shapes.Wall(
            normal=normal, dist=0.5 * self.agrid)
        self.lbf.add_boundary_from_shape(wall_shape)
        normal[dir_walls] = -1
        wall_shape = espressomd.shapes.Wall(
            normal=normal, dist=-(self.system.box_l[dir_walls] - 0.5 * self.agrid))
        self.lbf.add_boundary_from_shape(wall_shape)

        espressomd.utils.handle_errors("setup")

    @utx.skipIfMissingFeatures("EXTERNAL_FORCES")
    def test_ab_single_step(self):
        self.set_lb()
        self.system.part.clear()

        # Random velocities
        self.lbf[:, :, :].velocity = np.random.random((*self.lbf.shape, 3))
        force = [1, -2, 3]
        # Test several particle positions
        for pos in [[3 * self.agrid, 2 * self.agrid, 1 * self.agrid], [0, 0, 0],
                    self.system.box_l * 0.49,
                    self.system.box_l,
                    self.system.box_l * 0.99]:
            p = self.system.part.add(
                pos=pos,
                ext_force=force,
                propagation=espressomd.propagation.Propagation.TRANS_LB_TRACER)

            coupling_pos = p.pos
            # Nodes to which forces will be interpolated
            lb_nodes = tests_common.get_lb_nodes_around_pos(
                coupling_pos, self.lbf)

            np.testing.assert_allclose(
                [n.last_applied_force for n in lb_nodes],
                np.zeros((len(lb_nodes), 3)))
            self.system.integrator.run(1)

            v_fluid = np.copy(
                self.lbf.get_interpolated_velocity(
                    pos=coupling_pos))

            # Check particle velocity
            np.testing.assert_allclose(np.copy(p.v), v_fluid)

            # particle position
            np.testing.assert_allclose(
                np.copy(p.pos),
                coupling_pos + v_fluid * self.system.time_step)

            # check transfer of particle force to fluid
            applied_forces = np.array([n.last_applied_force for n in lb_nodes])
            np.testing.assert_allclose(
                np.sum(applied_forces, axis=0), force, atol=1E-10)

            # Check that last_applied_force gets cleared
            p.remove()
            self.system.integrator.run(1)
            applied_forces = np.array([n.last_applied_force for n in lb_nodes])
            np.testing.assert_allclose(
                np.sum(applied_forces, axis=0), [0, 0, 0])

    def test_advection(self):
        for direction in [0, 1, 2]:
            # System setup
            system = self.system

            # LB setup with walls
            ext_force = [0., 0., 0.]
            ext_force[direction] = 0.1
            dir_walls = (direction + 2) % 3
            box_l = 3 * [self.box_lw]
            box_l[dir_walls] = self.box_height
            system.box_l = box_l
            self.set_lb(ext_force_density=ext_force, dir_walls=dir_walls)

            # Establish steady state flow field
            system.integrator.run(400)

            # Add tracer in the fluid domain
            pos_initial = 3 * [3.5 * self.agrid]
            pos_initial[direction] = 0.5 * self.agrid
            p = system.part.add(
                pos=pos_initial,
                propagation=espressomd.propagation.Propagation.TRANS_LB_TRACER)

            # Perform integration
            system.time = 0
            for _ in range(2):
                system.integrator.run(100)
                # compute expected position
                lb_vel = self.lbf.get_interpolated_velocity(pos=p.pos)
                ref_dist = lb_vel[direction] * system.time
                tracer_dist = p.pos[direction] - pos_initial[direction]
                self.assertAlmostEqual(tracer_dist / ref_dist, 1., delta=0.01)

            system.lb = None
            system.part.clear()

    def test_zz_exceptions_without_lb(self):
        """
        Check behaviour without LB.
        """
        self.set_lb()
        system = self.system
        lbf = system.lb
        system.lb = None
        system.part.clear()
        p = system.part.add(pos=(0, 0, 0))
        with self.assertRaisesRegex(Exception, "The LB thermostat requires a LB fluid"):
            system.integrator.run(1)
        p.propagation = espressomd.propagation.Propagation.TRANS_LB_TRACER
        with self.assertRaisesRegex(Exception, "LB needs to be active for inertialess tracers"):
            system.integrator.run(1)
        system.lb = lbf
        self.system.thermostat.turn_off()
        with self.assertRaisesRegex(Exception, "The LB integrator requires the LB thermostat"):
            system.integrator.run(1)
