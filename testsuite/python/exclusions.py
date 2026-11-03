#
# Copyright (C) 2013-2026 The ESPResSo project
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
import unittest_decorators as utx
import espressomd
import espressomd.electrostatics

system = espressomd.System(box_l=[1.0, 1.0, 1.0])


@utx.skipIfMissingFeatures(['EXCLUSIONS'])
class Exclusions(ut.TestCase):
    system = system

    def setUp(self):
        self.system.box_l = 3 * [10]
        self.system.cell_system.skin = 0.4
        self.system.time_step = 0.01

    def tearDown(self):
        if espressomd.has_features("P3M"):
            self.system.electrostatics.clear()
        if espressomd.has_features("LENNARD_JONES"):
            self.system.non_bonded_inter[0, 0].lennard_jones.deactivate()
        self.system.part.clear()

    def test_add_remove(self):
        p0 = self.system.part.add(id=0, pos=[0, 0, 0])
        self.system.part.add(id=1, pos=[0, 0, 0])
        self.system.part.add(id=2, pos=[0, 0, 0])

        p0.add_exclusion(1)
        p0.add_exclusion(2)
        self.assertEqual(list(p0.exclusions), [1, 2])
        p0.delete_exclusion(1)
        self.assertEqual(list(p0.exclusions), [2])
        p0.delete_exclusion(2)
        self.assertEqual(list(p0.exclusions), [])

    def test_transfer(self):
        p0 = self.system.part.add(id=0, pos=[0, 0, 0], v=[1., 1., 1])
        self.system.part.add(id=1, pos=[0, 0, 0])
        self.system.part.add(id=2, pos=[0, 0, 0])
        self.system.part.add(id=3, pos=[0, 0, 0])

        p0.exclusions = [1, 2, 3]

        for _ in range(15):
            self.system.integrator.run(100)
            self.assertEqual(list(p0.exclusions), [1, 2, 3])

    @utx.skipIfMissingFeatures(['LENNARD_JONES'])
    def test_particle_property(self):
        self.system.non_bonded_inter[0, 0].lennard_jones.set_params(
            epsilon=1., sigma=2., cutoff=1.5, shift=0.0)

        p0 = self.system.part.add(id=0, pos=[0, 0, 0], type=0)
        p1 = self.system.part.add(id=1, pos=[1, 0, 0], type=0)

        pair_energy = self.system.analysis.energy()['total']
        self.assertGreater(pair_energy, 0.)

        pair_pressure = self.system.analysis.pressure()['total']
        self.assertGreater(pair_pressure, 0.)

        self.system.integrator.run(0)
        pair_force = p0.f[0]
        self.assertGreater(abs(pair_force), 0.)
        self.assertAlmostEqual(p1.f[0], -pair_force, places=7)

        p2 = self.system.part.add(id=2, pos=[2, 0, 0], type=0)
        self.system.integrator.run(0)
        self.assertAlmostEqual(self.system.analysis.energy()['total'],
                               2 * pair_energy)
        self.assertAlmostEqual(self.system.analysis.pressure()['total'],
                               2 * pair_pressure)
        self.assertAlmostEqual(p2.f[0], -pair_force, places=7)

        p1.exclusions = [0, 2]
        self.system.integrator.run(0)
        self.assertAlmostEqual(self.system.analysis.energy()['total'], 0)
        self.assertAlmostEqual(self.system.analysis.pressure()['total'], 0)
        self.assertAlmostEqual(p0.f[0], 0, places=7)
        self.assertAlmostEqual(p1.f[0], 0, places=7)
        self.assertAlmostEqual(p2.f[0], 0, places=7)

        p1.exclusions = [0]
        self.assertAlmostEqual(
            self.system.analysis.energy()['total'],
            pair_energy)
        self.assertAlmostEqual(self.system.analysis.pressure()['total'],
                               pair_pressure)
        self.system.integrator.run(0)
        self.assertAlmostEqual(p0.f[0], 0, places=7)
        self.assertAlmostEqual(p1.f[0], pair_force, places=7)
        self.assertAlmostEqual(p2.f[0], -pair_force, places=7)

        p1.exclusions = []
        self.assertAlmostEqual(self.system.analysis.energy()['total'],
                               2 * pair_energy)
        self.assertAlmostEqual(self.system.analysis.pressure()['total'],
                               2 * pair_pressure)
        self.system.integrator.run(0)
        self.assertAlmostEqual(p0.f[0], pair_force, places=7)
        self.assertAlmostEqual(p1.f[0], 0, places=7)
        self.assertAlmostEqual(p2.f[0], -pair_force, places=7)

        p1.exclusions = [0]
        self.assertAlmostEqual(
            self.system.analysis.energy()['total'],
            pair_energy)
        self.assertAlmostEqual(self.system.analysis.pressure()['total'],
                               pair_pressure)
        self.system.integrator.run(0)
        self.assertAlmostEqual(p0.f[0], 0, places=7)
        self.assertAlmostEqual(p1.f[0], pair_force, places=7)
        self.assertAlmostEqual(p2.f[0], -pair_force, places=7)

    @utx.skipIfMissingFeatures(['P3M'])
    def test_electrostatics_not_excluded(self):
        p0 = self.system.part.add(id=0, pos=[0, 0, 0], type=0, q=+1.)
        p1 = self.system.part.add(id=1, pos=[1, 0, 0], type=0, q=-1.)

        # Small alpha means large short-range contribution
        p3m = espressomd.electrostatics.P3M(
            prefactor=1, r_cut=3.0, accuracy=1e-3, mesh=32, cao=7, alpha=0.1,
            tune=False)
        self.system.electrostatics.solver = p3m

        # Only short-range part of the coulomb energy
        pair_energy = self.system.analysis.energy()[('coulomb', 0)]
        self.assertGreater(abs(pair_energy), 0.)

        self.system.integrator.run(0)
        pair_force = p0.f[0]
        self.assertGreater(abs(pair_force), 0.)
        self.assertAlmostEqual(p1.f[0], -pair_force, places=7)

        pair_pressure = self.system.analysis.pressure()[('coulomb', 0)]
        self.assertGreater(abs(pair_pressure), 0.)

        p0.exclusions = [1]
        # Force and energy should not be changed by the exclusion
        self.system.integrator.run(0)
        self.assertAlmostEqual(p0.f[0], pair_force, places=7)
        self.assertAlmostEqual(p1.f[0], -pair_force, places=7)
        self.assertAlmostEqual(self.system.analysis.energy()[('coulomb', 0)],
                               pair_energy, places=7)
        self.assertAlmostEqual(self.system.analysis.pressure()[('coulomb', 0)],
                               pair_pressure, places=7)


@utx.skipIfMissingFeatures(["EXCLUSIONS", "LENNARD_JONES"])
class ExclusionsGhostPairTest(ut.TestCase):
    system = system
    n_nodes = system.cell_system.get_state()["n_nodes"]
    # place the pair at the LJ minimum so the spurious contribution is
    # exactly -epsilon = -1.0 (at r = sigma it would be 0 and invisible)
    lj_min = 2.0**(1.0 / 6.0)

    def setUp(self):
        self.system.box_l = 3 * [10]
        self.system.cell_system.skin = 0.4
        self.system.time_step = 0.01
        self.system.non_bonded_inter[0, 0].lennard_jones.set_params(
            epsilon=1.0, sigma=1.0, cutoff=2.5, shift=0.0)
        self.node_grid = self.system.cell_system.node_grid

    def tearDown(self):
        self.system.part.clear()
        self.system.non_bonded_inter[0, 0].lennard_jones.deactivate()
        self.system.cell_system.node_grid = self.node_grid

    @ut.skipIf(n_nodes == 1, "only runs for 2+ MPI ranks (needs ghost layer)")
    def test_excluded_pair_cross_rank_ghost(self):
        """
        Excluded pair on different MPI ranks, interacting only through a
        cross-rank ghost image.  The ghost is rebuilt after the
        exclusion was set, so it carries an empty exclusion list.
        Since the pair is excluded, its contribution to the particle
        short-range energy must be zero; on buggy HEAD ``do_nonbonded``
        returns True for the (local, ghost) pair and the test fails
        with -1.0.
        """
        # split the box along x so the pair straddles the rank boundary
        self.system.cell_system.node_grid = [self.n_nodes, 1, 1]

        box_l = self.system.box_l[0]
        # minimum-image distance across the periodic x boundary
        # = 0.5 + (lj_min - 0.5) = lj_min
        pos_pair = ([0.5, 5.0, 5.0],
                    [box_l - (self.lj_min - 0.5), 5.0, 5.0])
        # far apart (min-image distance ~7 > cutoff): no interaction and
        # no ghost of one particle in the other's neighborhood
        pos_apart = ([2.5, 2.5, 5.0], [7.5, 7.5, 5.0])

        p0 = self.system.part.add(pos=pos_pair[0], type=0)
        p1 = self.system.part.add(pos=pos_pair[1], type=0)
        self.system.integrator.run(0)

        # sanity check: without exclusion the pair interacts through the
        # boundary with E = -epsilon = -1
        e_pair = self.system.analysis.particle_non_bonded_energy(p0)
        self.assertAlmostEqual(e_pair, -1.0, delta=1e-10,
                               msg="setup broken: pair not interacting")

        # Register the exclusion while the particles are far apart, then
        # move them back: the resort rebuilds the ghost copies through
        # the ghost communicator, which does not transfer exclusion
        # lists.  (Adding the exclusion while the ghost already exists
        # would patch the ghost in place and hide the bug.)
        p0.pos, p1.pos = pos_apart
        self.system.integrator.run(0)
        p0.add_exclusion(p1.id)
        p0.pos, p1.pos = pos_pair
        self.system.integrator.run(0)

        # correct behavior: the excluded pair contributes nothing.
        # On buggy HEAD do_nonbonded(local, ghost) returns True because
        # the ghost carries an empty exclusion list, yielding -1.0.
        for p in (p0, p1):
            energy = self.system.analysis.particle_non_bonded_energy(p)
            self.assertAlmostEqual(
                energy, 0.0, delta=1e-10,
                msg=f"excluded pair contributes {energy} to the short-range "
                f"energy of particle {p.id} (expected 0): the exclusion "
                "is ignored because the partner is a ghost with an empty "
                "exclusion list (OR instead of AND in do_nonbonded())")

    def test_excluded_local_pair_control(self):
        """
        Control: for a fully local excluded pair (both real particles in
        the same cell neighborhood, no ghost involved) both one-sided
        exclusion lists are populated, OR == AND, and the exclusion
        works.  This passes even on buggy HEAD and pins the failure of
        the ghost test on the OR combinator.
        """
        # separate along y so that under MPI (node grid split along x)
        # both particles stay on the same rank, far from any boundary
        p0 = self.system.part.add(pos=[2.5, 4.0, 5.0], type=0)
        p1 = self.system.part.add(pos=[2.5, 4.0 + self.lj_min, 5.0], type=0)
        p0.add_exclusion(p1.id)
        self.system.integrator.run(0)
        energy = self.system.analysis.particle_non_bonded_energy(p0)
        self.assertAlmostEqual(energy, 0.0, delta=1e-10)


if __name__ == "__main__":
    ut.main()
