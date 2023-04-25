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
import espressomd
import espressomd.shapes
import espressomd.lbboundaries
import espressomd.virtual_sites
import espressomd.utils


class VirtualSitesTracersCommon:
    box_height = 10.
    box_lw = 8.
    system = espressomd.System(box_l=(box_lw, box_lw, box_height))
    system.time_step = 0.05
    system.cell_system.skin = 0.1

    def setUp(self):
        self.system.box_l = (self.box_lw, self.box_lw, self.box_height)

    def tearDown(self):
        self.system.thermostat.turn_off()
        self.system.lbboundaries.clear()
        self.system.actors.clear()
        self.system.part.clear()

    def reset_lb(self, ext_force_density=(0, 0, 0), dir_walls=2):
        self.system.lbboundaries.clear()
        self.system.actors.clear()
        self.lbf = self.LBClass(
            kT=0.0, agrid=1, dens=1, visc=1.8,
            tau=self.system.time_step, ext_force_density=ext_force_density)
        self.system.actors.add(self.lbf)
        self.system.thermostat.set_lb(
            LB_fluid=self.lbf,
            act_on_virtual=False,
            gamma=1)

        # Setup boundaries
        normal = [0, 0, 0]
        normal[dir_walls] = 1
        walls = [espressomd.lbboundaries.LBBoundary() for k in range(2)]
        walls[0].set_params(shape=espressomd.shapes.Wall(
            normal=normal, dist=0.5))
        normal[dir_walls] = -1
        walls[1].set_params(shape=espressomd.shapes.Wall(
            normal=normal, dist=-(self.system.box_l[dir_walls] - 0.5)))

        for wall in walls:
            self.system.lbboundaries.add(wall)

        espressomd.utils.handle_errors("setup")

    def test_aa_method_switching(self):
        # Virtual sites should be disabled by default
        self.assertIsInstance(
            self.system.virtual_sites,
            espressomd.virtual_sites.VirtualSitesOff)

        # Switch implementation
        self.system.virtual_sites = espressomd.virtual_sites.VirtualSitesInertialessTracers()
        self.assertIsInstance(
            self.system.virtual_sites, espressomd.virtual_sites.VirtualSitesInertialessTracers)

    def test_advection(self):
        for direction in [0, 1, 2]:
            # System setup
            system = self.system
            system.virtual_sites = espressomd.virtual_sites.VirtualSitesInertialessTracers()

            # LB setup with walls
            ext_force = [0., 0., 0.]
            ext_force[direction] = 0.1
            dir_walls = (direction + 2) % 3
            box_l = 3 * [self.box_lw]
            box_l[dir_walls] = self.box_height
            system.box_l = box_l
            self.reset_lb(ext_force_density=ext_force, dir_walls=dir_walls)

            # Establish steady state flow field
            system.integrator.run(400)

            # Add tracer in the fluid domain
            pos_initial = [5.5, 5.5, 5.5]
            p = system.part.add(pos=pos_initial, virtual=True)

            # Perform integration
            system.time = 0
            for _ in range(2):
                system.integrator.run(100)
                # compute expected position
                lb_vel = self.lbf.get_interpolated_velocity(p.pos)
                ref_dist = lb_vel[direction] * system.time
                cur_dist = p.pos[direction] - pos_initial[direction]
                self.assertAlmostEqual(cur_dist / ref_dist, 1., delta=0.01)

            self.tearDown()

    def test_zz_without_lb(self):
        """Check behaviour without lb. Ignore non-virtual particles, complain on
        virtual ones.

        """
        self.reset_lb()
        system = self.system
        system.virtual_sites = espressomd.virtual_sites.VirtualSitesInertialessTracers()
        system.actors.clear()
        system.part.clear()
        p = system.part.add(pos=(0, 0, 0))
        system.integrator.run(1)
        p.virtual = True
        with self.assertRaises(Exception):
            system.integrator.run(1)
