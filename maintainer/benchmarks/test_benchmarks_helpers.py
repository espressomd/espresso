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

"""Pure-Python unit tests for benchmarks.py helpers (no espressomd needed)."""

import argparse
import pathlib
import tempfile
import numpy as np
import unittest as ut

import benchmarks


class StubCellSystem:
    def __init__(self, state):
        self._state = dict(state)

    def get_state(self):
        return dict(self._state)

    @property
    def node_grid(self):
        return self._state["node_grid"]

    @node_grid.setter
    def node_grid(self, value):
        self._state["node_grid"] = value

    @property
    def skin(self):
        return self._state.get("skin")

    @skin.setter
    def skin(self, value):
        self._state["skin"] = value


class StubSystem:
    def __init__(self, n_nodes=1, omp_num_threads=1, node_grid=None):
        self.cell_system = StubCellSystem({
            "n_nodes": n_nodes,
            "omp_num_threads": omp_num_threads,
            "node_grid": node_grid or [n_nodes, 1, 1],
        })


def make_args(**kwargs):
    ns = argparse.Namespace()
    defaults = dict(mode=None, state_file=None, skin=None,
                    particles_per_core=None, n_particles=None, _default_ppc=1000)
    defaults.update(kwargs)
    for k, v in defaults.items():
        setattr(ns, k, v)
    return ns


class Test(ut.TestCase):

    def test_n_cores_multiplies_ranks_and_threads(self):
        system = StubSystem(n_nodes=4, omp_num_threads=2)
        self.assertEqual(benchmarks.n_cores(system), 8)

    def test_resolve_n_part_particles_per_core_default(self):
        system = StubSystem(n_nodes=2, omp_num_threads=3)
        args = make_args(_default_ppc=100)
        self.assertEqual(benchmarks.resolve_n_part(system, args), 100 * 6)

    def test_resolve_n_part_explicit_ppc(self):
        system = StubSystem(n_nodes=2, omp_num_threads=1)
        args = make_args(particles_per_core=50)
        self.assertEqual(benchmarks.resolve_n_part(system, args), 100)

    def test_resolve_n_part_fixed_total(self):
        system = StubSystem(n_nodes=4, omp_num_threads=4)
        args = make_args(n_particles=777)
        self.assertEqual(benchmarks.resolve_n_part(system, args), 777)

    def test_validate_mode_requires_state_file(self):
        with self.assertRaises(SystemExit):
            benchmarks.validate_mode(make_args(mode="run", state_file=None))
        with self.assertRaises(SystemExit):
            benchmarks.validate_mode(make_args(mode="tune", state_file=None))
        benchmarks.validate_mode(make_args(mode=None))  # no exceptions

    def test_save_and_load_state_roundtrip(self):
        meta = {"skin": 0.4, "box_l": [10.0, 10.0, 10.0], "n_part": 3}
        pos = np.arange(9, dtype=float).reshape(3, 3)
        with tempfile.TemporaryDirectory() as d:
            path = str(pathlib.Path(d) / "state.npz")
            benchmarks.save_state(path, meta, pos=pos)
            meta_out, handle = benchmarks.load_state(path)
            self.assertEqual(meta_out, meta)
            np.testing.assert_array_equal(handle["pos"], pos)

    def test_topology_meta_and_verify_ok(self):
        system = StubSystem(n_nodes=2, omp_num_threads=2, node_grid=[2, 1, 1])
        meta = benchmarks.topology_meta(system)
        self.assertEqual(
            meta, {"n_nodes": 2, "node_grid": [2, 1, 1], "omp_num_threads": 2})
        # check against a fresh system with different node_grid but same counts
        other = StubSystem(n_nodes=2, omp_num_threads=2, node_grid=[1, 2, 1])
        benchmarks.verify_topology(other, meta)
        self.assertEqual(list(other.cell_system.node_grid), [2, 1, 1])

    def test_verify_topology_thread_mismatch(self):
        meta = {"n_nodes": 1, "node_grid": [1, 1, 1], "omp_num_threads": 4}
        system = StubSystem(n_nodes=1, omp_num_threads=2)
        with self.assertRaises(SystemExit):
            benchmarks.verify_topology(system, meta)

    def test_verify_topology_rank_mismatch(self):
        meta = {"n_nodes": 4, "node_grid": [4, 1, 1], "omp_num_threads": 1}
        system = StubSystem(n_nodes=2, omp_num_threads=1)
        with self.assertRaises(SystemExit):
            benchmarks.verify_topology(system, meta)

    def test_tune_skin_unless_fixed_uses_fixed_skin(self):
        system = StubSystem()
        args = make_args(skin=0.33)
        result = benchmarks.tune_skin_unless_fixed(system, args, 0.2, 1.0)
        self.assertAlmostEqual(result, 0.33, delta=1e-6)
        self.assertAlmostEqual(system.cell_system.skin, 0.33, delta=1e-6)

    def test_save_and_load_state_appends_npz_suffix(self):
        meta = {"n": 1}
        pos = np.zeros((2, 3))
        with tempfile.TemporaryDirectory() as d:
            stem = str(pathlib.Path(d) / "state")  # no .npz
            benchmarks.save_state(stem, meta, pos=pos)
            self.assertTrue(pathlib.Path(stem + ".npz").exists())
            meta_out, handle = benchmarks.load_state(stem)  # also no .npz
            self.assertEqual(meta_out, meta)
            np.testing.assert_array_equal(handle["pos"], pos)


if __name__ == "__main__":
    ut.main()
