#
# Copyright (C) 2026 The ESPResSo project
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
# MERCHANTABILITY or FITNESS FOR A PARTICLE PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
import unittest as ut
import unittest_decorators as utx
import espressomd
import espressomd.collision_detection
import espressomd.interactions
import numpy as np


@utx.skipIfMissingFeatures("COLLISION_DETECTION")
class BindCentersAcrossRanks(ut.TestCase):

    """
    Bond bookkeeping of the ``BindCenters`` collision protocol for a pair
    of particles owned by two different MPI ranks. The assertions state
    the behavior a bond created through ``system.bonds[bond].add()`` has;
    a failure means ``BindCenters`` deviates from it.
    """

    system = espressomd.System(box_l=[1.0, 1.0, 1.0])
    system.time_step = 0.001
    system.cell_system.skin = 0.05
    system.min_global_cut = 0.112
    # weak bond at its rest length: the pair stays put during the test
    bond = espressomd.interactions.HarmonicBond(k=1., r_0=0.02)
    system.bonded_inter.add(bond)

    default_node_grid = list(system.cell_system.node_grid)
    n_nodes = int(np.prod(default_node_grid))

    def setUp(self):
        if self.n_nodes < 2:
            self.skipTest("needs at least 2 MPI ranks")
        system = self.system
        system.collision_detection.protocol = espressomd.collision_detection.Off()
        system.part.clear()
        # Split the domain along z, so that the two particles straddling
        # the first internal boundary are owned by different ranks.
        system.cell_system.node_grid = [1, 1, self.n_nodes]
        boundary_z = system.box_l[2] / self.n_nodes
        self.p0 = system.part.add(pos=(0.5, 0.5, boundary_z - 0.01), id=0)
        self.p1 = system.part.add(pos=(0.5, 0.5, boundary_z + 0.01), id=1)
        self.assertNotEqual(self.p0.node, self.p1.node)

    def tearDown(self):
        system = self.system
        system.collision_detection.protocol = espressomd.collision_detection.Off()
        system.part.clear()
        system.cell_system.node_grid = self.default_node_grid

    def bind_by_collision(self):
        system = self.system
        system.collision_detection.protocol = espressomd.collision_detection.BindCenters(
            distance=0.05, bond_centers=self.bond)
        system.integrator.run(1, recalc_forces=True)
        system.collision_detection.protocol = espressomd.collision_detection.Off()

    def bind_by_interface(self):
        self.system.bonds[self.bond].add((self.p0, self.p1))

    def n_bonds(self):
        return len(self.p0.bonds) + len(self.p1.bonds)

    def refresh_ghosts(self):
        """Resort, so that ghost bond lists are overwritten by the bond
        lists of the real particles."""
        self.system.cell_system.resort()
        self.system.integrator.run(0, recalc_forces=True)

    def check_partner_removal(self, bind):
        """
        Removing the particle that does *not* own the bond must also
        remove the bond from its owner. This relies on the mirror entry
        stored on the real copy of the removed particle.
        """
        bind()
        self.assertGreaterEqual(self.n_bonds(), 1)
        self.refresh_ghosts()
        owner, partner = self.p0, self.p1
        if not owner.bonds:
            owner, partner = partner, owner
        if partner.bonds:
            self.skipTest("both particles own a bond (duplicate bond), "
                          "the mirror entry cannot be probed")
        partner.remove()
        self.assertEqual(owner.bonds, (), msg="dangling bond on the owner: "
                         "no mirror entry on the removed partner")
        # a dangling bond cannot be resolved by the force calculation
        self.system.integrator.run(1, recalc_forces=True)

    def test_partner_removal_reference(self):
        """Control: bond created through the bond interface."""
        self.check_partner_removal(self.bind_by_interface)

    def test_partner_removal_bind_centers(self):
        """Regression test: the mirror entry used to be written only on
        the ghost copy of the partner, never on its owning rank."""
        self.check_partner_removal(self.bind_by_collision)

    def test_no_duplicate_bond(self):
        """The two ranks must not both create the bond."""
        self.bind_by_collision()
        self.assertEqual(self.n_bonds(), 1)

    def test_no_duplicate_bond_later_steps(self):
        """The pair must not be bound again in later steps, neither
        before nor after the ghost bond lists were refreshed."""
        system = self.system
        system.collision_detection.protocol = espressomd.collision_detection.BindCenters(
            distance=0.05, bond_centers=self.bond)
        system.integrator.run(1, recalc_forces=True)
        self.assertEqual(self.n_bonds(), 1)
        system.integrator.run(5)
        self.assertEqual(self.n_bonds(), 1)
        self.refresh_ghosts()
        system.integrator.run(5)
        self.assertEqual(self.n_bonds(), 1)


if __name__ == "__main__":
    ut.main()
