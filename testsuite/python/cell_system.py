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
import numpy as np
import tests_common


class CellSystem(ut.TestCase):
    system = espressomd.System(box_l=[5., 5., 5.])
    n_nodes = system.cell_system.get_state()['n_nodes']

    def setUp(self):
        self.system.box_l = [5., 5., 5.]
        self.system.cell_system.skin = 0.

    def tearDown(self):
        self.system.part.clear()

    def test_cell_system(self):
        parameters = {
            "n_square": {"use_verlet_lists": False},
            "regular_decomposition": {"use_verlet_lists": True},
            "hybrid_decomposition": {"use_verlet_lists": False,
                                     "n_square_types": [1, 3, 5],
                                     "cutoff_regular": 1.27},
        }
        for cell_system, params_in in parameters.items():
            setter = getattr(self.system.cell_system, f"set_{cell_system}")
            setter(**params_in)
            params_in["decomposition_type"] = cell_system
            params_out = self.system.cell_system.get_params()
            tests_common.assert_params_match(self, params_in, params_out)
            params_out = self.system.cell_system.get_state()
            tests_common.assert_params_match(self, params_in, params_out)
        state = self.system.cell_system.get_state()
        self.assertIn("omp_num_threads", state)
        self.assertIsInstance(state["omp_num_threads"], int)
        self.assertGreaterEqual(state["omp_num_threads"], 1)

    def test_interface(self):
        self.system.cell_system.set_regular_decomposition()
        for value in [True, False]:
            self.system.cell_system.use_verlet_lists = value
            self.assertEqual(self.system.cell_system.use_verlet_lists, value)
        for value in [0.1, 0.]:
            self.system.cell_system.skin = value
            self.assertEqual(self.system.cell_system.skin, value)
        for key in ["decomposition_type", "n_square_types", "cutoff_regular",
                    "max_cut_nonbonded", "max_cut_bonded", "interaction_range"]:
            with self.assertRaisesRegex(RuntimeError, f"Parameter '{key}' is read-only"):
                setattr(self.system.cell_system, key, None)

    @ut.skipIf(n_nodes == 1, "Skipping test: only runs for n_nodes >= 2")
    def check_node_grid(self):
        system = self.system
        for i in range(3):
            node_grid_ref = [1, 1, 1]
            node_grid_ref[i] = self.n_nodes
            system.cell_system.node_grid = node_grid_ref
            node_grid = system.cell_system.get_state()['node_grid']
            np.testing.assert_array_equal(node_grid, node_grid_ref)

    def test_exceptions(self):
        system = self.system
        system.cell_system.skin = 0.1
        with self.assertRaisesRegex(ValueError, "Parameter 'skin' must be >= 0"):
            system.cell_system.skin = -2.
        self.assertAlmostEqual(system.cell_system.skin, 0.1, delta=1e-12)

        node_grid = system.cell_system.node_grid
        with self.assertRaisesRegex(RuntimeError, "Provided argument of type .+ is not convertible to 'Utils::Vector<int, 3>'"):
            system.cell_system.node_grid = [1, 2, 3, 4]
        np.testing.assert_array_equal(np.copy(system.cell_system.node_grid),
                                      np.copy(node_grid))
        with self.assertRaisesRegex(ValueError, rf"MPI world size {self.n_nodes} incompatible with new node grid \[1, 2, {self.n_nodes}\]"):
            system.cell_system.node_grid = [1, 2, self.n_nodes]
        np.testing.assert_array_equal(np.copy(system.cell_system.node_grid),
                                      np.copy(node_grid))

    def test_node_grid_regular(self):
        self.system.cell_system.set_regular_decomposition()
        self.check_node_grid()

    def test_node_grid_hybrid(self):
        self.system.cell_system.set_hybrid_decomposition(
            n_square_types={1}, cutoff_regular=0)
        self.check_node_grid()

    @utx.skipIfMissingFeatures(["LENNARD_JONES"])
    def test_force_survives_migration(self):
        """
        Regression test: a particle that moves to a new cell (or subdomain on
        a multi-rank run) during a global resort must carry its force with it.
        We compute a non-trivial pair force, move a particle far enough to
        leave its cell, then trigger a global resort WITHOUT recomputing forces
        (cell_system.resort()). The particle's force must be exactly the
        pre-resort value. Exercised for all three decompositions.
        """
        system = self.system
        system.part.clear()
        system.time_step = 0.01
        system.cell_system.skin = 0.1
        box_l = system.box_l[0]
        # cutoff (1.5) + skin (0.1) = 1.6 <= half box (2.5), so a regular
        # decomposition on one rank has a valid cell grid.
        system.non_bonded_inter[0, 0].lennard_jones.set_params(
            epsilon=1.0, sigma=1.0, cutoff=1.5, shift="auto")

        decompositions = [
            ("regular_decomposition", {}),
            ("n_square", {}),
            ("hybrid_decomposition",
             {"n_square_types": set(), "cutoff_regular": 1.5}),
        ]
        for name, kwargs in decompositions:
            with self.subTest(decomposition=name):
                system.part.clear()
                getattr(system.cell_system, f"set_{name}")(**kwargs)
                # A close-interacting pair so both feel a non-trivial force.
                p0 = system.part.add(pos=[0.5 * box_l, 0.5 * box_l,
                                          0.5 * box_l], type=0)
                p1 = system.part.add(pos=[0.5 * box_l + 1.1, 0.5 * box_l,
                                          0.5 * box_l], type=0)
                system.integrator.run(0)
                f0_before = np.copy(p0.f)
                f1_before = np.copy(p1.f)
                self.assertGreater(np.linalg.norm(f0_before), 1e-6)

                # Move p1 far across the box (and thus across any subdomain
                # boundary on a multi-rank run) so a global resort moves it.
                p1.pos = np.array([0.5 * box_l + 1.1, 0.5 * box_l,
                                   0.5 * box_l]) + np.array([0., 0., box_l / 2.])
                # Global resort moves the particle but does NOT recompute
                # forces. The force each particle carries must therefore be
                # exactly its pre-resort value.
                system.cell_system.resort()
                np.testing.assert_allclose(np.copy(p0.f), f0_before, atol=1e-12,
                                           err_msg=f"{name}: p0 force lost")
                np.testing.assert_allclose(np.copy(p1.f), f1_before, atol=1e-12,
                                           err_msg=f"{name}: p1 force lost")

        system.part.clear()
        system.non_bonded_inter[0, 0].lennard_jones.set_params(
            epsilon=0., sigma=0., cutoff=0., shift=0.)

    def test_rescale_isotropic_cubic_box(self):
        system = self.system
        system.box_l = [8., 8., 8.]
        p0, p1 = system.part.add(pos=[[1., 2., 48.], [-1., -2., 0.]])
        system.change_volume_and_rescale_particles(16., dir="xyz")
        np.testing.assert_almost_equal(np.copy(p0.pos), [2., 4., 96.])
        np.testing.assert_almost_equal(np.copy(p1.pos), [-2., -4., 0.])
        np.testing.assert_almost_equal(np.copy(system.box_l), [16., 16., 16.])
        system.change_volume_and_rescale_particles(8., dir="xyz")
        np.testing.assert_almost_equal(np.copy(p0.pos), [1., 2., 48.])
        np.testing.assert_almost_equal(np.copy(p1.pos), [-1., -2., 0.])
        np.testing.assert_almost_equal(np.copy(system.box_l), [8., 8., 8.])

    def test_rescale_isotropic_noncubic_box(self):
        system = self.system
        system.box_l = [4., 8., 12.]
        p0, p1 = system.part.add(pos=[[1., 2., 48.], [-48., -2., 0.]])
        system.change_volume_and_rescale_particles(8., dir="xyz")
        np.testing.assert_almost_equal(np.copy(p0.pos), [2., 2., 32.])
        np.testing.assert_almost_equal(np.copy(p1.pos), [-96., -2., 0.])
        np.testing.assert_almost_equal(np.copy(system.box_l), [8., 8., 8.])
        system.change_volume_and_rescale_particles(4., dir="xyz")
        np.testing.assert_almost_equal(np.copy(p0.pos), [1., 1., 16.])
        np.testing.assert_almost_equal(np.copy(p1.pos), [-48., -1., 0.])
        np.testing.assert_almost_equal(np.copy(system.box_l), [4., 4., 4.])

    def test_rescale_anisotropic(self):
        system = self.system
        system.box_l = [4., 4., 4.]
        p0, p1 = system.part.add(pos=[[1., 2., 3.], [-1., -2., 0.]])
        ref_pos = np.copy(system.part.all().pos)
        ref_box = np.array([4., 4., 4.])
        for i in range(3):
            system.change_volume_and_rescale_particles(8., dir=i)
            ref_pos[:, i] *= 2.
            ref_box[i] *= 2.
            np.testing.assert_almost_equal(np.copy(p0.pos), ref_pos[0])
            np.testing.assert_almost_equal(np.copy(p1.pos), ref_pos[1])
            np.testing.assert_almost_equal(np.copy(system.box_l), ref_box)
        for i in range(3):
            system.change_volume_and_rescale_particles(4., dir="xyz"[i])
            ref_pos[:, i] /= 2.
            ref_box[i] /= 2.
            np.testing.assert_almost_equal(np.copy(p0.pos), ref_pos[0])
            np.testing.assert_almost_equal(np.copy(p1.pos), ref_pos[1])
            np.testing.assert_almost_equal(np.copy(system.box_l), ref_box)

    @utx.skipIfMissingFeatures(["WCA"])
    @ut.skipIf(espressomd.has_features("FPE"),
               "cannot run with FPE instrumentation")
    def test_verlet_list_overflow(self):
        system = self.system
        # place all particles on top of each other
        system.part.add(pos=[[0, 0, 0]] * 1000)
        system.non_bonded_inter[0, 0].wca.set_params(epsilon=1., sigma=0.01)
        system.integrator.set_vv()

        system.cell_system.use_verlet_lists = True
        system.time_step = 0.01
        self.system.integrator.run(0)
        self.assertFalse(system.cell_system.use_verlet_lists)


if __name__ == "__main__":
    ut.main()
