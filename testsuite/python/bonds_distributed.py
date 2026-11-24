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
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
import unittest as ut
import espressomd
import espressomd.interactions
import numpy as np


class DistributedBonds(ut.TestCase):

    """
    Bonds are stored in the participating particles and move with them
    between MPI ranks. A periodic zigzag chain spans the box along x, which
    is split into one domain per rank, so that bonds cross every domain
    boundary and the periodic boundary.
    """

    system = espressomd.System(box_l=[12., 6., 6.])
    system.time_step = 0.01
    system.cell_system.skin = 0.2
    system.min_global_cut = 1.5
    n_nodes = system.cell_system.get_state()["n_nodes"]
    system.cell_system.node_grid = [n_nodes, 1, 1]
    harmonic = espressomd.interactions.HarmonicBond(k=10., r_0=0.5)
    angle = espressomd.interactions.AngleHarmonic(bend=2., phi0=np.pi)
    system.bonded_inter.add(harmonic)
    system.bonded_inter.add(angle)

    n_part = 24

    def setUp(self):
        x = 0.25 + 0.5 * np.arange(self.n_part)
        y = 3. + 0.2 * (-1.)**np.arange(self.n_part)
        pos = np.column_stack((x, y, np.full(self.n_part, 3.)))
        self.system.part.add(pos=pos)

    def tearDown(self):
        self.system.part.clear()
        self.system.cell_system.set_regular_decomposition()
        self.system.cell_system.node_grid = [self.n_nodes, 1, 1]

    def pair_rows(self):
        n = self.n_part
        return np.array([(i, (i + 1) % n) for i in range(n)])

    def angle_rows(self):
        n = self.n_part
        return np.array([(i, (i - 1) % n, (i + 1) % n) for i in range(n)])

    def add_chain_bonds(self):
        self.system.bonds[self.harmonic].add(self.pair_rows())
        self.system.bonds[self.angle].add(self.angle_rows())

    def reference_energy(self, pairs, angles):
        """Bonded energy, from the particle positions."""
        box = self.system.box_l
        pos = {p.id: p.pos for p in self.system.part.all()}

        def vec(i, j):
            d = pos[j] - pos[i]
            return d - box * np.round(d / box)

        energy = 0.
        for i, j in pairs:
            r = np.linalg.norm(vec(i, j))
            energy += 0.5 * self.harmonic.k * (r - self.harmonic.r_0)**2
        for vertex, a, b in angles:
            u, v = vec(vertex, a), vec(vertex, b)
            cos = np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v))
            phi = np.arccos(np.clip(cos, -1., 1.))
            energy += 0.5 * self.angle.bend * (phi - self.angle.phi0)**2
        return energy

    def owned_bonds(self):
        """All bonds, as listed on their owners."""
        return sorted((b[0]._bond_id, p.id, *b[1:])
                      for p in self.system.part.all() for b in p.bonds)

    def check_energy(self, pairs, angles):
        energy = self.system.analysis.energy()["bonded"]
        self.assertAlmostEqual(energy, self.reference_energy(pairs, angles),
                               delta=1e-10)

    def test_bonds_across_ranks(self):
        self.add_chain_bonds()
        self.check_energy(self.pair_rows(), self.angle_rows())
        self.assertEqual(len(self.owned_bonds()), 2 * self.n_part)
        # a bond is listed on its owner only
        p = self.system.part.by_id(5)
        self.assertEqual(p.bonds, ((self.harmonic, 6), (self.angle, 4, 6)))
        # the forces are consistent with the energy: they sum up to zero
        self.system.integrator.run(0)
        np.testing.assert_allclose(
            np.sum(self.system.part.all().f, axis=0), 0., atol=1e-10)

    def test_migration(self):
        self.add_chain_bonds()
        bonds = self.owned_bonds()
        energy = self.system.analysis.energy()["bonded"]
        # move every particle half a box along x, i.e. to another rank
        for shift in (6., 6.1, -3.):
            self.system.part.all().pos = self.system.part.all().pos + \
                [shift, 0., 0.]
            self.system.integrator.run(0)
            self.assertEqual(self.owned_bonds(), bonds)
            self.assertAlmostEqual(self.system.analysis.energy()["bonded"],
                                   energy, delta=1e-9)
        # migration during an integration
        self.system.part.all().v = [2., 0., 0.]
        self.system.integrator.run(400)
        self.assertEqual(self.owned_bonds(), bonds)
        self.check_energy(self.pair_rows(), self.angle_rows())

    def check_after_shuffle(self, bonds):
        """Move every particle to a random position, i.e. possibly several
        domains away, and back, with a global resort each time. The bonds
        are not evaluated in between, since the bonded particles are then
        too far apart."""
        rng = np.random.default_rng(42)
        pos = self.system.part.all().pos
        for _ in range(3):
            self.system.part.all().pos = rng.uniform(
                0., 1., (self.n_part, 3)) * self.system.box_l
            self.system.cell_system.resort()
            self.assertEqual(self.owned_bonds(), bonds)
        self.system.part.all().pos = pos
        self.system.cell_system.resort()
        self.system.integrator.run(0)
        self.assertEqual(self.owned_bonds(), bonds)
        self.check_energy(self.pair_rows(), self.angle_rows())

    def test_global_resort(self):
        self.add_chain_bonds()
        self.check_after_shuffle(self.owned_bonds())

    @ut.skipIf(n_nodes != 4, "requires 4 MPI ranks")
    def test_global_resort_2d_grid(self):
        # particles moving diagonally hop over an intermediate rank
        self.system.cell_system.node_grid = [2, 2, 1]
        self.add_chain_bonds()
        bonds = self.owned_bonds()
        self.check_after_shuffle(bonds)
        self.system.part.all().v = [2., 3., 0.]
        self.system.integrator.run(200)
        self.assertEqual(self.owned_bonds(), bonds)
        self.check_energy(self.pair_rows(), self.angle_rows())

    def test_other_decompositions(self):
        self.add_chain_bonds()
        bonds = self.owned_bonds()
        # switching the decomposition keeps the bonds
        self.system.cell_system.set_n_square()
        self.system.integrator.run(0)
        self.assertEqual(self.owned_bonds(), bonds)
        self.check_energy(self.pair_rows(), self.angle_rows())
        self.check_after_shuffle(bonds)
        self.system.part.all().v = [2., 1., 0.]
        self.system.integrator.run(200)
        self.assertEqual(self.owned_bonds(), bonds)
        self.check_energy(self.pair_rows(), self.angle_rows())
        # hybrid decomposition, with the chain in the N-square part, then in
        # the regular part (bonds between both parts cannot be resolved)
        self.system.part.all().type = 1
        self.system.cell_system.set_hybrid_decomposition(
            n_square_types={1}, cutoff_regular=1.)
        for chain_type in (1, 0):
            self.system.part.all().type = chain_type
            self.system.integrator.run(0)
            self.assertEqual(self.owned_bonds(), bonds)
            self.check_energy(self.pair_rows(), self.angle_rows())
            self.check_after_shuffle(bonds)
            self.system.integrator.run(200)
            self.assertEqual(self.owned_bonds(), bonds)
            self.check_energy(self.pair_rows(), self.angle_rows())
        # changes of the bonds after the particles have migrated
        self.system.bonds[self.harmonic].remove((3, 4))
        self.system.part.by_id(10).remove()
        pairs = [r for r in self.pair_rows()
                 if tuple(r) != (3, 4) and 10 not in r]
        angles = [r for r in self.angle_rows() if 10 not in r]
        self.check_energy(pairs, angles)

    def test_particle_removal(self):
        self.add_chain_bonds()
        # particle 5 owns bonds and is a partner in bonds owned by others,
        # which may reside on other ranks
        self.system.part.by_id(5).remove()
        pairs = [r for r in self.pair_rows() if 5 not in r]
        angles = [r for r in self.angle_rows() if 5 not in r]
        self.assertEqual(len(self.owned_bonds()), len(pairs) + len(angles))
        self.check_energy(pairs, angles)
        self.system.integrator.run(10)
        # a new particle with the same id does not inherit stale copies
        p = self.system.part.add(id=5, pos=[0.75, 3., 3.])
        self.assertEqual(p.bonds, ())
        self.system.bonds[self.harmonic].add((4, 5))
        self.system.bonds[self.harmonic].add((5, 6))
        self.assertEqual(self.system.part.by_id(4).bonds[0],
                         (self.harmonic, 5))
        self.system.integrator.run(10)

    def test_bulk_equals_single(self):
        self.add_chain_bonds()
        bulk_bonds = self.owned_bonds()
        bulk_energy = self.system.analysis.energy()["bonded"]
        self.system.part.clear()
        self.setUp()
        for row in self.pair_rows():
            self.system.bonds[self.harmonic].add(tuple(row))
        for row in self.angle_rows():
            self.system.bonds[self.angle].add(tuple(row))
        self.assertEqual(self.owned_bonds(), bulk_bonds)
        self.assertAlmostEqual(self.system.analysis.energy()["bonded"],
                               bulk_energy, delta=1e-10)
        # batches: a rejected batch leaves the previous batches in place
        self.system.part.clear()
        self.setUp()
        bonds = self.system.bonds[self.harmonic]
        batch_size = type(bonds)._batch_size
        try:
            type(bonds)._batch_size = 5
            rows = self.pair_rows()
            # 25 bonds in batches of 5: the duplicate is in the last batch
            with self.assertRaisesRegex(RuntimeError, "already exists"):
                bonds.add(np.vstack((rows, rows[:1])))
            self.assertEqual(len(self.owned_bonds()), 20)
        finally:
            type(bonds)._batch_size = batch_size

    def test_slice_setter_and_removal_across_ranks(self):
        self.add_chain_bonds()
        # replace the bonds owned by particles on all ranks at once
        everyone = self.system.part.all()
        everyone.bonds = [[(self.harmonic, (p.id + 2) % self.n_part)]
                          for p in everyone]
        pairs = [(i, (i + 2) % self.n_part) for i in range(self.n_part)]
        self.assertEqual(len(self.owned_bonds()), self.n_part)
        self.check_energy(pairs, [])
        # remove single bonds whose partners reside on other ranks
        for i in range(0, self.n_part, 3):
            self.system.bonds[self.harmonic].remove((i, (i + 2) % self.n_part))
        pairs = [r for r in pairs if r[0] % 3 != 0]
        self.check_energy(pairs, [])
        self.system.integrator.run(10)


if __name__ == "__main__":
    ut.main()
