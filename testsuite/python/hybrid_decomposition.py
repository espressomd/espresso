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
import unittest as ut
import unittest_decorators as utx
import tests_common
import espressomd
import itertools
import numpy as np


class HybridDecomposition(ut.TestCase):
    system = espressomd.System(box_l=3 * [50.0])
    original_node_grid = tuple(system.cell_system.node_grid)

    def setUp(self):
        self.system.cell_system.set_hybrid_decomposition(
            use_verlet_lists=False, n_square_types={1}, cutoff_regular=2.5)
        self.system.cell_system.node_grid = self.original_node_grid
        self.system.min_global_cut = 0.4
        self.system.time_step = 1e-3

    def tearDown(self):
        self.system.cell_system.set_regular_decomposition()
        self.system.part.clear()
        if espressomd.has_features(["ELECTROSTATICS"]):
            self.system.electrostatics.clear()
        if espressomd.has_features(["DIPOLES"]):
            self.system.magnetostatics.clear()

    def check_resort(self):
        n_part = 2352

        # number of particles that should end up in the respective child
        # decomposition
        n_n_square = n_part // 4
        n_regular = n_part - n_n_square

        # Add the particles on node 0, so that they have to be resorted
        particles = self.system.part.add(
            pos=n_part * [(0, 0, 0)], type=n_n_square * [0, 0, 0, 1])
        parts_per_decomposition = self.system.cell_system.get_state()[
            'parts_per_decomposition']
        self.assertEqual(parts_per_decomposition['n_square'], n_n_square)
        self.assertEqual(parts_per_decomposition['regular'], n_regular)

        # And now change their positions
        particles.pos = self.system.box_l * \
            np.random.random((n_part, 3))

        # Add an interacting particle in a corner of the box
        self.system.part.add(pos=(0.01, 0.01, 0.01), type=0)
        if espressomd.has_features(['LENNARD_JONES']):
            self.system.non_bonded_inter[0, 1].lennard_jones.set_params(
                epsilon=1.0, sigma=3.0, cutoff=6.0, shift=0.1)
            ref_energy = self.system.analysis.energy()['total']
            assert ref_energy > 10.

        # Distribute the particles on the nodes
        part_dist = self.system.cell_system.resort()

        # Check that we did not lose particles
        self.assertEqual(sum(part_dist), n_part + 1)

        # Check that we can still access all the particles
        # This basically checks if part_node and local_particles
        # is still in a valid state after the particle exchange
        self.assertEqual(len(self.system.part.all()), n_part + 1)
        self.assertEqual(sum(self.system.part.all().type), n_n_square)

        # Check that the system is still valid
        if espressomd.has_features(['LENNARD_JONES']):
            # energy calculation
            new_energy = self.system.analysis.energy()['total']
            self.assertEqual(new_energy, ref_energy)
        # force calculation
        self.system.integrator.run(0, recalc_forces=True)

    def prepare_hybrid_setup(self, n_part_small=0, n_part_large=0):
        """Setup system with small and large particles, minimize
         energy and setup thermostat. Particles have random
         initial velocities; cutoff of small particles is 2.5.

        """
        assert espressomd.has_features(["LENNARD_JONES"])
        box_l = self.system.box_l[0]
        if n_part_small > 0:
            self.system.part.add(
                type=[0] * n_part_small,
                pos=np.random.random((n_part_small, 3)) * box_l,
                v=np.random.randn(n_part_small, 3))
        if n_part_large > 0:
            self.system.part.add(
                type=[1] * n_part_large,
                pos=np.random.random((n_part_large, 3)) * box_l,
                v=np.random.randn(n_part_large, 3))
        self.assertEqual(len(self.system.part.all()),
                         n_part_small + n_part_large)

        # setup interactions
        self.system.non_bonded_inter[0, 0].lennard_jones.set_params(
            epsilon=1, sigma=1, cutoff=2.5, shift="auto")
        # mixing rule: sigma = 0.5 * (s_large + s_small)
        self.system.non_bonded_inter[0, 1].lennard_jones.set_params(
            epsilon=1, sigma=2.5, cutoff=2**(1 / 6) * 2.5, shift="auto")
        self.system.non_bonded_inter[1, 1].lennard_jones.set_params(
            epsilon=1, sigma=4, cutoff=10, shift="auto")

        # remove overlap
        self.system.integrator.set_steepest_descent(
            f_max=0, gamma=30, max_displacement=0.01)
        self.system.integrator.run(0)
        old_force = np.max(np.linalg.norm(self.system.part.all().f, axis=-1))
        while old_force > n_part_small + n_part_large:
            self.system.integrator.run(20)
            force = np.max(np.linalg.norm(self.system.part.all().f, axis=-1))
            old_force = force

        self.system.integrator.set_vv()

    def test_resort(self):
        self.check_resort()

    @ut.skipIf(system.cell_system.get_state()["n_nodes"] != 4,
               "Skipping test: only runs for n_nodes >= 4")
    def test_resort_alternating(self):
        # check particle resorting when the left and right cells are different
        self.system.cell_system.node_grid = [4, 1, 1]
        self.check_resort()

    def test_position_rounding(self):
        """This places a particle on the box boundary,
           with parameters that could cause problems with
           rounding."""
        self.system.box_l = [50.0, 50.0, 50.0]
        self.system.cell_system.skin = 0.4
        self.system.min_global_cut = 12.0 / 4.25
        self.system.part.add(pos=[25, 25, 0])
        self.assertEqual(1, len(self.system.part))

    @utx.skipIfMissingFeatures(["LENNARD_JONES"])
    def test_non_bonded_loop_trace(self):
        self.prepare_hybrid_setup(n_part_small=50, n_part_large=50)
        cutoff = 2.5
        tests_common.check_non_bonded_loop_trace(self, self.system, cutoff)

    @utx.skipIfMissingFeatures(["LENNARD_JONES"])
    def test_against_nsquare(self):
        self.prepare_hybrid_setup(n_part_small=150, n_part_large=50)

        self.run_comparison(self.system.cell_system.set_n_square)

    def run_comparison(self, set_reference_decomposition_fun):
        steps_per_round = 20
        for _ in range(40):
            # integrate using hybrid and calculate energy and forces
            self.system.cell_system.set_hybrid_decomposition(
                n_square_types={1}, cutoff_regular=2.5)
            self.system.integrator.run(steps_per_round)
            energy_mix = self.system.analysis.energy()
            forces_mix = np.copy(self.system.part.all().f)

            # compare to reference cell system for consistency
            set_reference_decomposition_fun()
            self.system.integrator.run(0)

            energy_ref = self.system.analysis.energy()
            forces_ref = np.copy(self.system.part.all().f)
            self.assertAlmostEqual(
                energy_mix["non_bonded"],
                energy_ref["non_bonded"],
                delta=1e-9)
            self.assertAlmostEqual(
                energy_mix[("non_bonded", 0, 0)],
                energy_ref[("non_bonded", 0, 0)],
                delta=1e-9)
            self.assertAlmostEqual(
                energy_mix[("non_bonded", 0, 1)],
                energy_ref[("non_bonded", 0, 1)],
                delta=1e-9)
            self.assertAlmostEqual(
                energy_mix[("non_bonded", 1, 1)],
                energy_ref[("non_bonded", 1, 1)],
                delta=1e-9)
            np.testing.assert_allclose(
                forces_mix, forces_ref, rtol=1e-9, atol=0.)

    def test_sort_into_child_decs(self):
        """Assert that particles end up in the respective child
        decomposition, depending on their type. Also, check that
        changing particle type or n_square_types results in the
        expected resort.

        """
        n_parts = 3
        parts = self.system.part.add(
            pos=np.random.random((n_parts, 3)) * self.system.box_l[0],
            type=np.random.randint(2, size=n_parts))
        for ndx, types in enumerate(itertools.product([0, 1], repeat=n_parts)):
            parts.type = types
            n_square_type = ndx % 2
            n_n_square = n_parts - \
                np.sum(types) if n_square_type == 0 else np.sum(types)
            n_regular = n_parts - n_n_square

            self.system.cell_system.set_hybrid_decomposition(
                n_square_types={n_square_type}, cutoff_regular=0)
            parts_per_decomposition = self.system.cell_system.get_state()[
                'parts_per_decomposition']
            self.assertEqual(parts_per_decomposition['n_square'], n_n_square)
            self.assertEqual(parts_per_decomposition['regular'], n_regular)

    def valid_p3m_parameters(self):
        return {"r_cut": 16.11914, "alpha": 0.14735, "mesh": 12, "cao": 5,
                "tune": False, "accuracy": 1e-5, "prefactor": 1.0}

    def valid_dp3m_parameters(self):
        return {"r_cut": 17.8418, "alpha": 0.09919, "mesh": 4, "cao": 2,
                "tune": False, "accuracy": 1e-5, "prefactor": 1.0}

    def add_particles(self, charge=False, dip=False):
        eps = 1e-14
        corner = np.full(3, eps)
        edge = [self.system.box_l[0] - eps,
                self.system.box_l[1] - eps,
                self.system.box_l[2] / 2]
        p = self.system.part.add(
            pos=[corner, [0.5, 2.6, 7.6], edge],
            type=[0, 1, 0],
        )
        if charge:
            p.q = np.array([-1., 1., 0.])
        if dip:
            p.dip = np.array([(1., 0., 0.), (0., -1., 0.), (0., 0., 1.)])

    @ut.skipIf(system.cell_system.get_state()["n_nodes"] != 1,
               "Skipping test: only runs for n_nodes == 1")
    @utx.skipIfMissingFeatures(["P3M", "LENNARD_JONES"])
    def test_against_regular_p3m(self):
        import espressomd.electrostatics

        self.prepare_hybrid_setup(n_part_small=0, n_part_large=0)
        self.add_particles(charge=True)
        coulomb_interaction = espressomd.electrostatics.P3M(
            **self.valid_p3m_parameters()
        )
        self.system.electrostatics.solver = coulomb_interaction

        self.run_comparison(self.system.cell_system.set_regular_decomposition)

    @ut.skipIf(system.cell_system.get_state()["n_nodes"] != 1,
               "Skipping test: only runs for n_nodes == 1")
    @utx.skipIfMissingFeatures(["DP3M", "LENNARD_JONES"])
    def test_against_regular_dp3m(self):
        import espressomd.magnetostatics

        self.prepare_hybrid_setup(n_part_small=0, n_part_large=0)
        self.add_particles(dip=True)
        actor = espressomd.magnetostatics.DipolarP3M(
            **self.valid_dp3m_parameters()
        )
        self.system.magnetostatics.solver = actor

        self.run_comparison(self.system.cell_system.set_regular_decomposition)

    @ut.skipIf(system.cell_system.get_state()["n_nodes"] in [1, 3],
               "Skipping test: only runs for n_nodes not in [1, 3]")
    @utx.skipIfMissingFeatures(["P3M"])
    def test_mpi_exception_p3m(self):
        import espressomd.electrostatics

        self.system.cell_system.set_regular_decomposition()
        self.system.analysis.energy()
        actor = espressomd.electrostatics.P3M(
            **self.valid_p3m_parameters()
        )
        self.system.electrostatics.solver = actor

        with self.assertRaises(Exception):
            self.system.cell_system.set_hybrid_decomposition(
                use_verlet_lists=False, n_square_types={1}, cutoff_regular=2.5)

    @ut.skipIf(system.cell_system.get_state()["n_nodes"] in [1, 3],
               "Skipping test: only runs for n_nodes not in [1, 3]")
    @utx.skipIfMissingFeatures(["DP3M"])
    def test_mpi_exception_dp3m(self):
        import espressomd.magnetostatics

        self.system.cell_system.set_regular_decomposition()
        self.system.analysis.energy()
        actor = espressomd.magnetostatics.DipolarP3M(
            **self.valid_dp3m_parameters()
        )
        self.system.magnetostatics.solver = actor

        with self.assertRaises(Exception):
            self.system.cell_system.set_hybrid_decomposition(
                use_verlet_lists=False, n_square_types={1}, cutoff_regular=2.5)


if __name__ == "__main__":
    ut.main()
