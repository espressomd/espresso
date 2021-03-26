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
import numpy as np

import espressomd
from espressomd import shapes, lbboundaries
import espressomd.interactions
try:
    from espressomd.virtual_sites import VirtualSitesInertialessTracers, VirtualSitesOff
except ImportError:
    pass
from espressomd.utils import handle_errors
from tests_common import get_lb_nodes_around_pos

import unittest_decorators as utx


class VirtualSitesTracersCommon:
    box_height = 10.
    box_lw = 8.
    system = espressomd.System(box_l=(box_lw, box_lw, box_height))
    system.time_step = 0.08
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
        walls = [lbboundaries.LBBoundary() for k in range(2)]
        walls[0].set_params(shape=shapes.Wall(normal=[0, 0, 1], dist=0.5))
        walls[1].set_params(shape=shapes.Wall(normal=[0, 0, -1],
                                              dist=-self.box_height - 0.5))

        for wall in walls:
            self.system.lbboundaries.add(wall)

        handle_errors("setup")

    def test_aa_method_switching(self):
        # Virtual sites should be disabled by default
        self.assertIsInstance(self.system.virtual_sites, VirtualSitesOff)

        # Switch implementation
        self.system.virtual_sites = VirtualSitesInertialessTracers()
        self.assertIsInstance(
            self.system.virtual_sites, VirtualSitesInertialessTracers)

    @utx.skipIfMissingFeatures("EXTERNAL_FORCES")
    def test_ab_single_step(self):
        self.reset_lb()
        self.system.lbboundaries.clear()
        self.system.part.clear()
        self.system.virtual_sites = VirtualSitesInertialessTracers()

        # Random velocities
        for n in self.lbf.nodes():
            n.velocity = np.random.random(3) - .5
        force = [1, -2, 3]
        # Test several particle positions
        for pos in [[3, 2, 1], [0, 0, 0],
                    self.system.box_l * 0.49,
                    self.system.box_l,
                    self.system.box_l * 0.99]:
            p = self.system.part.add(pos=pos, ext_force=force, virtual=True)

            coupling_pos = p.pos
            # Nodes to which forces will be interpolated
            lb_nodes = get_lb_nodes_around_pos(
                coupling_pos, self.lbf)

            np.testing.assert_allclose(
                [n.last_applied_force for n in lb_nodes],
                np.zeros((len(lb_nodes), 3)))
            self.system.integrator.run(1)

            v_fluid = np.copy(self.lbf.get_interpolated_velocity(coupling_pos))

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
        self.reset_lb(ext_force_density=[0.1, 0, 0])
        # System setup
        system = self.system

        system.virtual_sites = VirtualSitesInertialessTracers()

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
            self.assertAlmostEqual(p.pos[0] / dist, 1, delta=0.001)

    def compute_dihedral_angle(self, pos0, pos1, pos2, pos3):
        # first normal vector
        n1 = np.cross((pos1 - pos0), (pos2 - pos0))
        n2 = np.cross((pos2 - pos0), (pos3 - pos0))

        norm1 = np.linalg.norm(n1)
        norm2 = np.linalg.norm(n2)
        n1 = n1 / norm1
        n2 = n2 / norm2

        cos_alpha = min(1, np.dot(n1, n2))
        alpha = np.arccos(cos_alpha)
        return alpha

    def test_tribend(self):
        # two triangles with bending interaction
        # move nodes, should relax back

        system = self.system
        system.virtual_sites = VirtualSitesInertialessTracers()
        system.thermostat.set_langevin(kT=0, gamma=10, seed=1)

        # Add four particles
        p0 = system.part.add(pos=[5, 5, 5])
        p1 = system.part.add(pos=[5, 5, 6])
        p2 = system.part.add(pos=[5, 6, 6])
        p3 = system.part.add(pos=[5, 6, 5])

        # Add first triel, weak modulus
        tri1 = espressomd.interactions.IBM_Triel(
            ind1=p0.id, ind2=p1.id, ind3=p2.id, elasticLaw="Skalak", k1=0.1, k2=0, maxDist=2.4)
        system.bonded_inter.add(tri1)
        p0.add_bond((tri1, p1, p2))

        # Add second triel, strong modulus
        tri2 = espressomd.interactions.IBM_Triel(
            ind1=p0.id, ind2=p2.id, ind3=p3.id, elasticLaw="Skalak", k1=10, k2=0, maxDist=2.4)
        system.bonded_inter.add(tri2)
        p0.add_bond((tri2, p2, p3))

        # Add bending
        tribend = espressomd.interactions.IBM_Tribend(
            ind1=p0.id, ind2=p1.id, ind3=p2.id, ind4=p3.id, kb=1, refShape="Initial")
        system.bonded_inter.add(tribend)
        p0.add_bond((tribend, p1, p2, p3))

        # twist
        system.part[:].pos = system.part[:].pos + np.random.random((4, 3))

        # Perform integration
        system.integrator.run(200)
        angle = self.compute_dihedral_angle(p0.pos, p1.pos, p2.pos, p3.pos)
        self.assertLess(angle, 1.5E-3)

    def test_triel(self):
        system = self.system
        system.virtual_sites = VirtualSitesInertialessTracers()
        system.thermostat.set_langevin(kT=0, gamma=1, seed=1)

        # Add particles: 0-2 are not bonded, 3-5 are bonded
        non_bound = system.part.add(pos=[[5, 5, 5], [5, 5, 6], [5, 6, 6]])

        p3 = system.part.add(pos=[2, 5, 5])
        p4 = system.part.add(pos=[2, 5, 6])
        p5 = system.part.add(pos=[2, 6, 6])

        # Add triel for 3-5
        tri = espressomd.interactions.IBM_Triel(
            ind1=p3.id, ind2=p4.id, ind3=p5.id, elasticLaw="Skalak", k1=15,
            k2=0, maxDist=2.4)
        system.bonded_inter.add(tri)
        p3.add_bond((tri, p4, p5))

        system.part[:].pos = system.part[:].pos + np.array((
            (0, 0, 0), (1, -.2, .3), (1, 1, 1),
            (0, 0, 0), (1, -.2, .3), (1, 1, 1)))

        distorted_pos = np.copy(non_bound.pos)

        system.integrator.run(110)
        dist1bound = system.distance(p3, p4)
        dist2bound = system.distance(p3, p5)

        # check bonded particles. Distance should restore to initial config
        self.assertAlmostEqual(dist1bound, 1, delta=0.02)
        self.assertAlmostEqual(dist2bound, np.sqrt(2), delta=0.01)

        # check not bonded particles. Positions should still be distorted
        np.testing.assert_allclose(np.copy(non_bound.pos), distorted_pos)

    def test_zz_without_lb(self):
        """Check behaviour without lb. Ignore non-virtual particles, complain on
        virtual ones.

        """
        self.reset_lb()
        system = self.system
        system.virtual_sites = VirtualSitesInertialessTracers()
        system.actors.clear()
        system.part.clear()
        p = system.part.add(pos=(0, 0, 0))
        system.integrator.run(1)
        p.virtual = True
        with self.assertRaises(Exception):
            system.integrator.run(1)
