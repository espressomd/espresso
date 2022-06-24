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

    def tearDown(self):
        self.system.part.clear()
        self.system.actors.clear()
        self.system.thermostat.turn_off()

    def reset_lb(self, ext_force_density=(0, 0, 0)):
        self.system.lbboundaries.clear()
        self.lbf = self.LBClass(
            kT=0.0, agrid=1, dens=1, visc=1.8,
            tau=self.system.time_step, ext_force_density=ext_force_density)
        self.system.actors.add(self.lbf)
        self.system.thermostat.set_lb(
            LB_fluid=self.lbf,
            act_on_virtual=False,
            gamma=1)

        # Setup boundaries
        walls = [espressomd.lbboundaries.LBBoundary() for k in range(2)]
        walls[0].set_params(shape=espressomd.shapes.Wall(
            normal=[0, 0, 1], dist=0.5))
        walls[1].set_params(shape=espressomd.shapes.Wall(
            normal=[0, 0, -1], dist=-self.box_height - 0.5))

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
        self.reset_lb(ext_force_density=[0.1, 0, 0])
        # System setup
        system = self.system

        system.virtual_sites = espressomd.virtual_sites.VirtualSitesInertialessTracers()

        # Establish steady state flow field
        p = system.part.add(pos=(0, 5.5, 5.5), virtual=True)
        system.integrator.run(400)

        p.pos = (0, 5.5, 5.5)
        system.time = 0

        # Perform integration
        for _ in range(2):
            system.integrator.run(100)
            # compute expected position
            dist = self.lbf.get_interpolated_velocity(p.pos)[0] * system.time
            self.assertAlmostEqual(p.pos[0] / dist, 1, delta=0.005)

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
