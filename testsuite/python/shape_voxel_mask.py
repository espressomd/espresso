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

import contextlib
import itertools
import math
import os
import pathlib
import pickle
import re
import tempfile
import unittest as ut
import unittest_decorators as utx
import numpy as np

import espressomd
import espressomd.lb
import espressomd.shapes

with contextlib.suppress(ImportError):
    import h5py  # h5py has to be imported *after* espressomd (MPI)


def node_grids(n_nodes):
    """All node grids with ``n_nodes`` ranks."""
    for grid in itertools.product(range(1, n_nodes + 1), repeat=3):
        if math.prod(grid) == n_nodes:
            yield grid


def solid_cubes(mask, voxel_size, origin, box_l):
    """
    Lower corners of the solid voxel cubes, including periodic images.
    Mask voxels whose center lies outside of the box are ignored.
    """
    indices = np.argwhere(mask)
    origin = np.asarray(origin, dtype=float)
    box_l = np.asarray(box_l, dtype=float)
    centers = origin + (indices + 0.5) * voxel_size
    in_box = np.all((centers >= 0.) & (centers < box_l), axis=1)
    lower = origin + indices[in_box] * voxel_size
    images = np.array(list(itertools.product((-1, 0, 1), repeat=3)))
    return (lower[np.newaxis, :, :] + (images * box_l)
            [:, np.newaxis, :]).reshape((-1, 3))


def reference_distance(pos, cubes, voxel_size):
    """Distance from ``pos`` to the nearest solid cube."""
    closest = np.clip(pos, cubes, cubes + voxel_size)
    return np.min(np.linalg.norm(pos - closest, axis=1))


def reference_is_inside(pos, mask, voxel_size, origin, box_l):
    """Solid state of the voxel that contains ``pos``."""
    index = []
    for axis in range(3):
        j = math.floor((pos[axis] - origin[axis]) / voxel_size)
        center = origin[axis] + (j + 0.5) * voxel_size
        center -= math.floor(center / box_l[axis]) * box_l[axis]
        m = math.floor((center - origin[axis]) / voxel_size)
        if m < 0 or m >= mask.shape[axis]:
            return False
        index.append(m)
    return bool(mask[tuple(index)])


@utx.skipIfMissingFeatures(["H5MD"])
@utx.skipIfMissingModules("h5py")
class VoxelMaskTest(ut.TestCase):
    system = espressomd.System(box_l=[16., 16., 16.])
    n_nodes = system.cell_system.get_state()["n_nodes"]

    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.node_grid = np.copy(self.system.cell_system.node_grid)

    def tearDown(self):
        self.system.cell_system.node_grid = self.node_grid
        self.temp_dir.cleanup()

    def write_mask(self, mask, name="mask.h5", dtype=np.uint8):
        path = pathlib.Path(self.temp_dir.name) / name
        with h5py.File(path, "w") as h5:
            h5.create_dataset("mask", data=np.asarray(mask, dtype=dtype))
        return path

    def check_distances(self, mask, voxel_size, origin, halo, n_samples):
        box_l = np.copy(self.system.box_l)
        path = self.write_mask(mask)
        shape = espressomd.shapes.VoxelMask(
            file=str(path), dataset="mask", voxel_size=voxel_size,
            origin=origin, halo=halo)
        cubes = solid_cubes(mask, voxel_size, origin, box_l)
        rng = np.random.default_rng(seed=42)
        n_pore = 0
        n_solid = 0
        while n_pore < n_samples:
            pos = rng.random(3) * box_l
            if reference_is_inside(pos, mask, voxel_size, origin, box_l):
                self.assertTrue(shape.is_inside(position=pos))
                n_solid += 1
                continue
            ref = reference_distance(pos, cubes, voxel_size)
            dist, vec = shape.calc_distance(position=pos)
            self.assertAlmostEqual(dist, ref, delta=1e-12)
            self.assertAlmostEqual(np.linalg.norm(vec), ref, delta=1e-12)
            self.assertFalse(shape.is_inside(position=pos))
            n_pore += 1
        self.assertGreater(n_solid, 0)

    def test_distance_matches_reference(self):
        self.system.box_l = [16., 16., 16.]
        rng = np.random.default_rng(seed=17)
        mask = rng.random((16, 16, 16)) < 0.3
        self.check_distances(mask, voxel_size=1., origin=(0., 0., 0.),
                             halo=3, n_samples=200)

    def test_distance_fine_voxels_shifted_origin(self):
        self.system.box_l = [16., 16., 16.]
        rng = np.random.default_rng(seed=23)
        mask = rng.random((32, 32, 32)) < 0.3
        self.check_distances(mask, voxel_size=0.5, origin=(-0.25, 0.5, 0.),
                             halo=4, n_samples=100)

    def test_orientation(self):
        self.system.box_l = [8., 12., 16.]
        mask = np.zeros((8, 12, 16), dtype=bool)
        mask[1, 5, 11] = True
        path = self.write_mask(mask, dtype=bool)
        # positional arguments
        shape = espressomd.shapes.VoxelMask(str(path), "mask", 1.)
        np.testing.assert_array_equal(shape.mask_shape, [8, 12, 16])
        self.assertTrue(shape.is_inside(position=[1.5, 5.5, 11.5]))
        self.assertFalse(shape.is_inside(position=[11.5, 5.5, 1.5]))
        self.assertFalse(shape.is_inside(position=[1.5, 11.5, 5.5]))

    def test_mask_smaller_than_box(self):
        self.system.box_l = [12., 12., 12.]
        mask = np.ones((8, 12, 12), dtype=np.uint8)
        path = self.write_mask(mask)
        shape = espressomd.shapes.VoxelMask(
            file=str(path), dataset="mask", voxel_size=1., halo=3)
        self.assertTrue(shape.is_inside(position=[4., 1., 1.]))
        self.assertFalse(shape.is_inside(position=[10., 1., 1.]))
        dist, _ = shape.calc_distance(position=[10., 1., 1.])
        self.assertAlmostEqual(dist, 2., delta=1e-12)
        dist, _ = shape.calc_distance(position=[10., 6., 6.])
        self.assertAlmostEqual(dist, 2., delta=1e-12)

    def test_mask_larger_than_box(self):
        self.system.box_l = [16., 16., 16.]
        mask = np.zeros((20, 16, 16), dtype=np.uint8)
        mask[18, :, :] = 1
        path = self.write_mask(mask)
        shape = espressomd.shapes.VoxelMask(
            file=str(path), dataset="mask", voxel_size=1.)
        for x in np.arange(0.5, 16., 1.):
            for pos in ([x, 0.5, 0.5], [x, 7.5, 12.5]):
                self.assertFalse(shape.is_inside(position=pos))
                dist, _ = shape.calc_distance(position=pos)
                self.assertGreater(dist, shape.halo * shape.voxel_size)

    def test_local_memory(self):
        box_l = 32.
        halo = 1
        self.system.box_l = 3 * [box_l]
        mask = np.zeros((32, 32, 32), dtype=np.uint8)
        mask[::3, ::5, ::7] = 1
        path = self.write_mask(mask)
        shape = espressomd.shapes.VoxelMask(
            file=str(path), dataset="mask", voxel_size=1., halo=halo)
        for node_grid in node_grids(self.n_nodes):
            self.system.cell_system.node_grid = node_grid
            expected = []
            for node_pos in itertools.product(*(range(n) for n in node_grid)):
                count = 1
                for axis in range(3):
                    local_l = box_l / node_grid[axis]
                    left = node_pos[axis] * local_l
                    right = (node_pos[axis] + 1) * local_l
                    lower = math.floor(left) - halo
                    upper = math.ceil(right) + halo
                    count *= upper - lower
                expected.append(count)
            counts = shape.get_local_voxel_counts()
            memory = shape.get_local_memory_bytes()
            self.assertEqual(len(counts), self.n_nodes)
            self.assertEqual(len(memory), self.n_nodes)
            self.assertEqual(sorted(counts), sorted(expected))
            for i in range(self.n_nodes):
                if self.n_nodes > 1:
                    self.assertLess(counts[i], mask.size)
                self.assertLessEqual(memory[i], counts[i] / 8 + 64)
            # queries work after the node grid change
            self.assertTrue(shape.is_inside(position=[0.5, 0.5, 0.5]))
            self.assertTrue(shape.is_inside(position=[30.5, 25.5, 28.5]))
            self.assertFalse(shape.is_inside(position=[31.5, 25.5, 28.5]))

    def test_rasterize_raises(self):
        mask = np.ones((16, 16, 16), dtype=np.uint8)
        path = self.write_mask(mask)
        shape = espressomd.shapes.VoxelMask(
            file=str(path), dataset="mask", voxel_size=1.)
        with self.assertRaisesRegex(RuntimeError, "cannot be rasterized"):
            shape.call_method("rasterize", grid_size=[16, 16, 16],
                              grid_spacing=1., grid_offset=0.5)

    def test_missing_file_raises(self):
        path = pathlib.Path(self.temp_dir.name) / "missing.h5"
        stderr_path = pathlib.Path(self.temp_dir.name) / "stderr.txt"
        # capture the C-level standard error of the head rank
        stderr_fd = os.dup(2)
        try:
            with open(stderr_path, "w") as stderr_file:
                os.dup2(stderr_file.fileno(), 2)
                message = re.escape(f"Cannot open HDF5 file '{path}': ")
                with self.assertRaisesRegex(RuntimeError, message + ".+"):
                    espressomd.shapes.VoxelMask(
                        file=str(path), dataset="mask", voxel_size=1.)
        finally:
            os.dup2(stderr_fd, 2)
            os.close(stderr_fd)
        self.assertNotIn("HDF5", stderr_path.read_text())

    def test_missing_dataset_raises(self):
        path = self.write_mask(np.ones((4, 4, 4)))
        with self.assertRaisesRegex(RuntimeError, "Cannot find dataset 'xy'"):
            espressomd.shapes.VoxelMask(
                file=str(path), dataset="xy", voxel_size=1.)

    def test_invalid_parameters_raise(self):
        path = self.write_mask(np.ones((16, 16, 16)))
        with self.assertRaisesRegex(ValueError, "'halo' must be >= 1"):
            espressomd.shapes.VoxelMask(
                file=str(path), dataset="mask", voxel_size=1., halo=0)
        with self.assertRaisesRegex(ValueError, "'voxel_size' must be > 0"):
            espressomd.shapes.VoxelMask(
                file=str(path), dataset="mask", voxel_size=0.)

    def test_union_raises(self):
        path = self.write_mask(np.ones((16, 16, 16)))
        shape = espressomd.shapes.VoxelMask(
            file=str(path), dataset="mask", voxel_size=1.)
        wall = espressomd.shapes.Wall(normal=[1., 0., 0.], dist=1.)
        union = espressomd.shapes.Union()
        union.add(wall)
        message = "A VoxelMask cannot be part of a Union"
        with self.assertRaisesRegex(ValueError, message):
            union.add(shape)
        with self.assertRaisesRegex(ValueError, message):
            union.add([shape])
        with self.assertRaisesRegex(ValueError, message):
            espressomd.shapes.Union(_objects=[wall, shape])
        # the union is still usable
        self.assertEqual(union.size(), 1)
        dist, _ = union.calc_distance(position=[3., 1., 1.])
        self.assertAlmostEqual(dist, 2., delta=1e-12)

    def test_parameters(self):
        path = self.write_mask(np.ones((16, 15, 14)))
        shape = espressomd.shapes.VoxelMask(
            file=str(path), dataset="mask", voxel_size=1.,
            origin=[0.5, 0., -0.5], halo=3)
        self.assertEqual(shape.file, str(path))
        self.assertEqual(shape.dataset, "mask")
        self.assertEqual(shape.voxel_size, 1.)
        np.testing.assert_array_equal(shape.origin, [0.5, 0., -0.5])
        self.assertEqual(shape.halo, 3)
        np.testing.assert_array_equal(shape.mask_shape, [16, 15, 14])
        with self.assertRaisesRegex(RuntimeError, "read-only"):
            shape.halo = 4

    def test_pickle(self):
        self.system.box_l = [16., 16., 16.]
        rng = np.random.default_rng(seed=5)
        mask = rng.random((16, 16, 16)) < 0.3
        path = self.write_mask(mask)
        shape = espressomd.shapes.VoxelMask(
            file=str(path), dataset="mask", voxel_size=1., halo=3)
        shape_copy = pickle.loads(pickle.dumps(shape))
        self.assertEqual(shape_copy.file, str(path))
        self.assertEqual(shape_copy.halo, 3)
        for pos in rng.random((10, 3)) * 16.:
            dist_ref, vec_ref = shape.calc_distance(position=pos)
            dist, vec = shape_copy.calc_distance(position=pos)
            self.assertEqual(dist, dist_ref)
            np.testing.assert_array_equal(np.copy(vec), np.copy(vec_ref))


def wca_force(dist, sigma=1., epsilon=1.):
    """Magnitude of the repulsive WCA force at distance ``dist``."""
    if dist >= 2.**(1. / 6.) * sigma:
        return 0.
    sr6 = (sigma / dist)**6
    return 24. * epsilon * (2. * sr6**2 - sr6) / dist


def wca_energy(dist, sigma=1., epsilon=1.):
    """Repulsive WCA energy at distance ``dist``."""
    if dist >= 2.**(1. / 6.) * sigma:
        return 0.
    sr6 = (sigma / dist)**6
    return 4. * epsilon * (sr6**2 - sr6) + epsilon


def reference_force(pos, cubes, voxel_size):
    """WCA force of the solid cubes on a particle at ``pos``."""
    closest = np.clip(pos, cubes, cubes + voxel_size)
    diff = pos - closest
    dists = np.linalg.norm(diff, axis=1)
    i = np.argmin(dists)
    return wca_force(dists[i]) * diff[i] / dists[i]


@utx.skipIfMissingFeatures(["H5MD", "LENNARD_JONES"])
@utx.skipIfMissingModules("h5py")
class VoxelMaskConstraint(ut.TestCase):
    system = VoxelMaskTest.system
    n_nodes = VoxelMaskTest.n_nodes

    def setUp(self):
        self.system.time_step = 0.01
        self.system.cell_system.skin = 0.4
        self.system.box_l = [16., 16., 16.]
        self.temp_dir = tempfile.TemporaryDirectory()
        self.node_grid = np.copy(self.system.cell_system.node_grid)
        self.system.non_bonded_inter[0, 1].wca.set_params(
            epsilon=1., sigma=1.)

    def tearDown(self):
        self.system.constraints.clear()
        self.system.part.clear()
        self.system.periodicity = [True, True, True]
        self.system.cell_system.node_grid = self.node_grid
        self.system.non_bonded_inter.reset()
        self.temp_dir.cleanup()

    def write_mask(self, mask, name="mask.h5"):
        path = pathlib.Path(self.temp_dir.name) / name
        with h5py.File(path, "w") as h5:
            h5.create_dataset("mask", data=np.asarray(mask, dtype=np.uint8))
        return path

    def add_constraint(self, mask, halo=3):
        path = self.write_mask(mask)
        shape = espressomd.shapes.VoxelMask(
            file=str(path), dataset="mask", voxel_size=1.,
            origin=(0., 0., 0.), halo=halo)
        self.system.constraints.add(
            shape=shape, particle_type=1, penetrable=False)
        return shape

    def setup_porous(self):
        rng = np.random.default_rng(seed=17)
        mask = rng.random((16, 16, 16)) < 0.3
        self.mask = mask
        self.add_constraint(mask)
        cubes = solid_cubes(mask, 1., (0., 0., 0.), [16., 16., 16.])
        positions = []
        forces = []
        while len(positions) < 100:
            pos = rng.random(3) * 16.
            if reference_is_inside(pos, mask, 1., (0., 0., 0.), [16.] * 3):
                continue
            if reference_distance(pos, cubes, 1.) <= 0.6:
                continue
            positions.append(pos)
            forces.append(reference_force(pos, cubes, 1.))
        self.system.part.add(pos=positions, type=[0] * len(positions))
        self.assertGreater(np.count_nonzero(np.linalg.norm(forces, axis=1)), 5)
        return np.array(positions), np.array(forces)

    def check_forces(self, forces):
        self.system.integrator.run(0, recalc_forces=True)
        # particle ids are consecutive and ordered like the positions
        np.testing.assert_allclose(
            np.copy(self.system.part.all().f), forces, rtol=0., atol=1e-10)

    def test_forces_match_reference(self):
        _, forces = self.setup_porous()
        self.check_forces(forces)

    def test_energy_matches_reference(self):
        positions, _ = self.setup_porous()
        cubes = solid_cubes(self.mask, 1., (0., 0., 0.), [16.] * 3)
        energy = sum(wca_energy(reference_distance(pos, cubes, 1.))
                     for pos in positions)
        self.assertGreater(energy, 0.)
        energies = self.system.analysis.energy()
        self.assertAlmostEqual(energies["non_bonded", 0, 1], energy,
                               delta=1e-10)
        self.assertAlmostEqual(energies["total"], energy, delta=1e-10)

    def test_forces_after_node_grid_change(self):
        # No shape method may be called after a node grid change: it would
        # reload the data and mask a missing lifecycle hook.
        _, forces = self.setup_porous()
        for grid in node_grids(self.n_nodes):
            self.system.cell_system.node_grid = grid
            self.check_forces(forces)

    def test_reload_error(self):
        mask = np.zeros((16, 16, 16), dtype=bool)
        mask[15, :, :] = True
        shape = self.add_constraint(mask)
        part = self.system.part.add(pos=[0.5, 8., 8.], type=0)
        path = pathlib.Path(shape.file)
        moved = path.with_name("moved.h5")
        message = re.escape(f"Cannot open HDF5 file '{path}'")
        path.rename(moved)
        # every rank fails to reload, and all ranks raise the error
        with self.assertRaisesRegex(RuntimeError, message):
            self.system.periodicity = [False, True, True]
        # the data is discarded, every use retries the reload
        with self.assertRaisesRegex(Exception, message):
            self.system.integrator.run(0, recalc_forces=True)
        with self.assertRaisesRegex(RuntimeError, message):
            shape.calc_distance(position=[0.5, 8., 8.])
        if self.n_nodes > 1:
            old_grid = np.copy(self.system.cell_system.node_grid)
            grid = [g for g in node_grids(self.n_nodes)
                    if not np.array_equal(g, old_grid)][0]
            with self.assertRaisesRegex(RuntimeError, message):
                self.system.cell_system.node_grid = grid
            np.testing.assert_array_equal(
                np.copy(self.system.cell_system.node_grid), old_grid)
        moved.rename(path)
        self.system.integrator.run(0, recalc_forces=True)
        np.testing.assert_array_equal(np.copy(part.f), 0.)
        self.system.periodicity = [True, True, True]
        self.system.integrator.run(0, recalc_forces=True)
        self.assertGreater(part.f[0], 1.)

    def test_periodicity_change(self):
        mask = np.zeros((16, 16, 16), dtype=bool)
        mask[15, :, :] = True
        self.add_constraint(mask)
        part = self.system.part.add(pos=[0.5, 8., 8.], type=0)
        self.system.integrator.run(0, recalc_forces=True)
        f_periodic = np.copy(part.f)
        self.assertGreater(f_periodic[0], 1.)
        np.testing.assert_allclose(f_periodic[1:], 0., atol=1e-10)
        self.system.periodicity = [False, True, True]
        self.system.integrator.run(0, recalc_forces=True)
        np.testing.assert_array_equal(np.copy(part.f), 0.)
        self.system.periodicity = [True, True, True]
        self.system.integrator.run(0, recalc_forces=True)
        np.testing.assert_allclose(np.copy(part.f), f_periodic, atol=1e-10)

    def make_shape(self, mask, halo=3):
        path = self.write_mask(mask)
        return espressomd.shapes.VoxelMask(
            file=str(path), dataset="mask", voxel_size=1.,
            origin=(0., 0., 0.), halo=halo)

    def test_add_after_periodicity_change(self):
        # the shape is not part of a constraint when the periodicity changes
        mask = np.zeros((16, 16, 16), dtype=bool)
        mask[15, :, :] = True
        shape = self.make_shape(mask)
        self.system.periodicity = [False, True, True]
        self.system.constraints.add(shape=shape, particle_type=1)
        part = self.system.part.add(pos=[0.5, 8., 8.], type=0)
        self.system.integrator.run(0, recalc_forces=True)
        np.testing.assert_array_equal(np.copy(part.f), 0.)

    def test_add_after_box_change(self):
        # the shape is not part of a constraint when the box changes,
        # e.g. a shape created before the system is set up
        mask = np.zeros((8, 8, 8), dtype=bool)
        mask[:2, :, :] = True
        shape = self.make_shape(mask)
        self.system.box_l = [8., 8., 8.]
        self.system.constraints.add(shape=shape, particle_type=1)
        p1 = self.system.part.add(pos=[3., 4., 4.], type=0)
        p2 = self.system.part.add(pos=[7., 4., 4.], type=0)
        self.system.integrator.run(0, recalc_forces=True)
        f_ref = wca_force(1.)
        np.testing.assert_allclose(np.copy(p1.f), [f_ref, 0., 0.], atol=1e-10)
        np.testing.assert_allclose(np.copy(p2.f), [-f_ref, 0., 0.], atol=1e-10)

    def test_far_from_interface(self):
        # the LJ cutoff reaches beyond the halo
        self.system.non_bonded_inter[0, 1].lennard_jones.set_params(
            epsilon=1., sigma=1., cutoff=2.5, shift=0.)
        mask = np.zeros((16, 16, 16), dtype=bool)
        mask[0, :, :] = True
        shape = self.make_shape(mask, halo=1)
        constraint = self.system.constraints.add(shape=shape, particle_type=1)
        part = self.system.part.add(pos=[8., 8., 8.], type=0)
        # no interface within the halo: no phantom interaction
        self.assertEqual(self.system.analysis.energy()["non_bonded", 0, 1], 0.)
        # the integrator rejects a halo smaller than max_cutoff + skin / 2
        message = (r"halo is too small .* particle type 1: .* got "
                   r"halo \* voxel_size = 1, max_cutoff \+ skin / 2 = 2\.7")
        with self.assertRaisesRegex(Exception, message):
            self.system.integrator.run(0, recalc_forces=True)
        # a sufficient halo works
        self.system.constraints.remove(constraint)
        shape = self.make_shape(mask, halo=3)
        self.system.constraints.add(shape=shape, particle_type=1)
        self.system.integrator.run(0, recalc_forces=True)
        np.testing.assert_array_equal(np.copy(part.f), 0.)
        self.assertEqual(self.system.analysis.energy()["non_bonded", 0, 1], 0.)
        part.pos = [2.5, 8., 8.]
        self.system.integrator.run(0, recalc_forces=True)
        self.assertLess(part.f[0], 0.)

    def test_force_zero_far_from_wall(self):
        mask = np.zeros((16, 16, 16), dtype=bool)
        mask[7:9, 7:9, 7:9] = True
        self.add_constraint(mask)
        part = self.system.part.add(pos=[2.5, 2.5, 2.5], type=0)
        self.system.integrator.run(0, recalc_forces=True)
        np.testing.assert_array_equal(np.copy(part.f), [0., 0., 0.])


@utx.skipIfMissingFeatures(["H5MD", "WALBERLA"])
@utx.skipIfMissingModules("h5py")
class VoxelMaskLB(ut.TestCase):
    system = VoxelMaskTest.system
    n_nodes = VoxelMaskTest.n_nodes

    def setUp(self):
        self.system.time_step = 0.01
        self.system.cell_system.skin = 0.4
        self.system.box_l = [16., 16., 16.]
        self.temp_dir = tempfile.TemporaryDirectory()
        self.node_grid = np.copy(self.system.cell_system.node_grid)

    def tearDown(self):
        self.system.lb = None
        self.system.box_l = [16., 16., 16.]
        self.system.cell_system.node_grid = self.node_grid
        self.temp_dir.cleanup()

    def write_mask(self, mask, name="mask.h5"):
        path = pathlib.Path(self.temp_dir.name) / name
        with h5py.File(path, "w") as h5:
            h5.create_dataset("mask", data=np.asarray(mask, dtype=np.uint8))
        return path

    def make_shape(self, mask, origin=(0., 0., 0.), halo=3):
        path = self.write_mask(mask)
        return espressomd.shapes.VoxelMask(
            file=str(path), dataset="mask", voxel_size=1., origin=origin,
            halo=halo)

    def make_lb(self, n_ghost_layers=1):
        lattice = espressomd.lb.Lattice(
            agrid=1., n_ghost_layers=n_ghost_layers)
        lbf = espressomd.lb.LBFluid(
            lattice=lattice, kinematic_viscosity=1., density=1., tau=0.01)
        self.system.lb = lbf
        return lbf

    def porous_mask(self, n_voxels):
        rng = np.random.default_rng(seed=17)
        return rng.random(3 * (n_voxels,)) < 0.3

    def compatible_node_grids(self, n_cells):
        """Node grids that split the LB lattice evenly."""
        return [grid for grid in node_grids(self.n_nodes)
                if all(n_cells % n == 0 for n in grid)]

    def check_boundary_flags(self, mask):
        shape = self.make_shape(mask)
        lbf = self.make_lb()
        lbf.add_boundary_from_shape(shape)
        np.testing.assert_array_equal(
            np.copy(lbf[:, :, :].is_boundary), mask.astype(bool))
        self.assertGreater(np.count_nonzero(mask), 0)

    def test_boundary_flags_match_mask(self):
        grids = self.compatible_node_grids(16)
        if not grids:
            self.skipTest("no node grid splits a 16^3 lattice evenly")
        self.system.cell_system.node_grid = grids[0]
        self.check_boundary_flags(self.porous_mask(16))

    def test_boundary_flags_all_node_grids(self):
        # 24 cells can be split evenly by every node grid of up to 4 ranks
        self.system.box_l = [24., 24., 24.]
        mask = self.porous_mask(24)
        grids = self.compatible_node_grids(24)
        self.assertGreater(len(grids), 0)
        for grid in grids:
            with self.subTest(node_grid=grid):
                self.system.cell_system.node_grid = grid
                try:
                    self.check_boundary_flags(mask)
                finally:
                    self.system.lb = None

    def test_alignment_warning(self):
        self.system.box_l = [24., 24., 24.]
        self.system.cell_system.node_grid = self.compatible_node_grids(24)[0]
        mask = np.zeros((24, 24, 24), dtype=bool)
        mask[4:8, 4:8, 4:8] = True
        shape = self.make_shape(mask, origin=(0.5, 0., 0.))
        lbf = self.make_lb()
        with self.assertWarns(UserWarning):
            lbf.add_boundary_from_shape(shape)

    def test_halo_smaller_than_ghost_layers_raises(self):
        self.system.box_l = [24., 24., 24.]
        self.system.cell_system.node_grid = self.compatible_node_grids(24)[0]
        mask = self.porous_mask(24)
        lbf = self.make_lb(n_ghost_layers=3)
        # the ghost-layer node centers reach 2.5 voxels beyond the local box
        shape = self.make_shape(mask, halo=2)
        message = (r"VoxelMask halo is too small for the LB ghost layers.*"
                   r"halo=2, voxel_size=1, agrid=1, n_ghost_layers=3")
        with self.assertRaisesRegex(ValueError, message):
            lbf.add_boundary_from_shape(shape)
        np.testing.assert_array_equal(np.copy(lbf[:, :, :].is_boundary), False)
        # a halo that covers the ghost layers works
        shape = self.make_shape(mask, halo=3)
        lbf.add_boundary_from_shape(shape)
        np.testing.assert_array_equal(
            np.copy(lbf[:, :, :].is_boundary), mask.astype(bool))

    def test_reload_error(self):
        mask = np.zeros((24, 24, 24), dtype=bool)
        mask[4:8, 4:8, 4:8] = True
        shape = self.make_shape(mask)
        path = pathlib.Path(shape.file)
        moved = path.with_name("moved.h5")
        path.rename(moved)
        # the box change makes the shape data stale
        self.system.box_l = [24., 24., 24.]
        self.system.cell_system.node_grid = self.compatible_node_grids(24)[0]
        lbf = self.make_lb()
        message = re.escape(f"Cannot open HDF5 file '{path}'")
        with self.assertRaisesRegex(RuntimeError, message):
            lbf.add_boundary_from_shape(shape)
        np.testing.assert_array_equal(np.copy(lbf[:, :, :].is_boundary), False)
        moved.rename(path)
        lbf.add_boundary_from_shape(shape)
        np.testing.assert_array_equal(
            np.copy(lbf[:, :, :].is_boundary), mask.astype(bool))

    def test_array_velocity_raises(self):
        self.system.box_l = [24., 24., 24.]
        self.system.cell_system.node_grid = self.compatible_node_grids(24)[0]
        mask = np.zeros((24, 24, 24), dtype=bool)
        mask[4:8, 4:8, 4:8] = True
        shape = self.make_shape(mask)
        lbf = self.make_lb()
        velocity = np.zeros(tuple(lbf.shape) + (3,))
        with self.assertRaisesRegex(ValueError, "VoxelMask"):
            lbf.add_boundary_from_shape(shape, velocity)
        # the boundary force needs a global grid of the shape
        with self.assertRaisesRegex(RuntimeError, "cannot be rasterized"):
            lbf.get_boundary_force_from_shape(shape)


if __name__ == "__main__":
    ut.main()
