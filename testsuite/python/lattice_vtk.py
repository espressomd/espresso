#
# Copyright (C) 2010-2026 The ESPResSo project
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

import os
import pathlib
import tempfile
import contextlib
import numpy as np

import espressomd
import espressomd.lb
import espressomd.shapes
import espressomd.electrokinetics

with contextlib.suppress(ImportError):
    import espressomd.io.vtk


class TestVTK:
    system = espressomd.System(box_l=[6, 7, 3])
    system.time_step = 0.1
    system.cell_system.skin = 0.4

    lattice_params = {}

    def setUp(self):
        self.lattice = self.lattice_class(
            n_ghost_layers=2, agrid=0.5, **self.lattice_params)
        self.actor = self.add_actor()

    def tearDown(self):
        self.clear_actors()

    @ut.skipIf(system.cell_system.get_state()["n_nodes"] > 4,
               "this test is slow on more than 4 MPI ranks")
    def test_exceptions(self):
        label_invalid_obs = f"test_vtk_{self.vtk_id}_invalid_obs"
        error_msg = rf"Only the following VTK observables are supported: \[{repr(sorted(self.valid_obs))[1:-1]}\], got 'dens'"  # nopep8
        with self.assertRaisesRegex(ValueError, error_msg):
            self.vtk_class(
                identifier=label_invalid_obs, delta_N=0, observables=["dens"])
        vtk_manual_id = f"test_vtk_{self.vtk_id}_manual"
        vtk_auto_id = f"test_vtk_{self.vtk_id}_auto"
        vtk_manual = self.vtk_class(
            identifier=vtk_manual_id, delta_N=0, observables=["density"])
        vtk_auto = self.vtk_class(
            identifier=vtk_auto_id, delta_N=1, observables=["density"])
        self.actor.add_vtk_writer(vtk=vtk_manual)
        self.actor.add_vtk_writer(vtk=vtk_auto)
        with self.assertRaisesRegex(RuntimeError, "Automatic VTK callbacks cannot be triggered manually"):
            vtk_auto.write()
        with self.assertRaisesRegex(RuntimeError, "Manual VTK callbacks cannot be disabled"):
            vtk_manual.disable()
        with self.assertRaisesRegex(RuntimeError, "Manual VTK callbacks cannot be enabled"):
            vtk_manual.enable()
        with self.assertRaisesRegex(RuntimeError, "already exists"):
            self.actor.add_vtk_writer(vtk=self.vtk_class(
                identifier=vtk_manual_id, delta_N=0, observables=[]))
        with self.assertRaisesRegex(RuntimeError, "already attached to this lattice"):
            self.actor.add_vtk_writer(vtk=self.actor.vtk_writers[0])
        with self.assertRaisesRegex(RuntimeError, "not attached to this lattice"):
            self.actor.remove_vtk_writer(vtk=self.vtk_class(
                identifier=vtk_manual_id, delta_N=0, observables=[]))
        with self.assertRaisesRegex(RuntimeError, "Cannot attach VTK object to multiple lattices"):
            self.make_actor().add_vtk_writer(vtk=vtk_manual)
        with self.assertRaisesRegex(RuntimeError, "Detached VTK objects cannot be attached again"):
            self.actor.remove_vtk_writer(vtk=vtk_manual)
            self.actor.add_vtk_writer(vtk=vtk_manual)
        with self.assertRaisesRegex(ValueError, "Parameter 'delta_N' must be >= 0"):
            self.vtk_class(identifier="a", delta_N=-1, observables=[])
        with self.assertRaisesRegex(ValueError, "Parameter 'identifier' cannot be empty"):
            self.vtk_class(identifier="", delta_N=0, observables=[])
        with self.assertRaisesRegex(ValueError, "cannot be a filepath"):
            self.vtk_class(
                identifier=f"test{os.sep}test", delta_N=0, observables=[])

        # can still use VTK when the actor has been cleared but not deleted
        label_cleared = f"test_vtk_{self.vtk_id}_cleared"
        vtk_cleared = self.vtk_class(
            identifier=label_cleared, observables=["density"])
        self.actor.add_vtk_writer(vtk=vtk_cleared)
        self.clear_actors()
        vtk_cleared.write()

        # cannot use VTK when no lattice is attached to it
        label_unattached = f"test_vtk_{self.vtk_id}_unattached"
        label_unattached = self.vtk_class(
            identifier=label_unattached, observables=[])
        with self.assertRaisesRegex(RuntimeError, "This VTK object isn't attached to a lattice"):
            label_unattached.write()

    @utx.skipIfMissingModules("espressomd.io.vtk")
    def test_exceptions_invalid_files(self):
        with tempfile.TemporaryDirectory() as tmp_directory:
            root = pathlib.Path(tmp_directory)
            invalid_vtk_file = root / "invalid_file.vtu"
            invalid_vtk_file.write_text(1000 * "\n    ")
            with self.assertRaisesRegex(RuntimeError, "is not a compliant XML file"):
                espressomd.io.vtk.VTKReader().parse(invalid_vtk_file)
            invalid_vtk_file.write_text('<VTKFile type="UnknownGrid"/>')
            with self.assertRaisesRegex(NotImplementedError, "Unknown VTK file format 'UnknownGrid'"):
                espressomd.io.vtk.VTKReader().parse(invalid_vtk_file)


class TestLBVTK(TestVTK):
    include_boundaries = True

    valid_obs = ["density", "velocity_vector", "pressure_tensor", "boundary"]

    def write_obs(self):
        obs = ["density", "velocity_vector", "pressure_tensor"]
        if self.include_boundaries:
            obs.append("boundary")
        return obs

    def make_actor(self):
        return self.lb_class(
            lattice=self.lattice, tau=0.1, density=1.2, kinematic_viscosity=1.,
            ext_force_density=[0., 0.03, 0.], **self.lb_params)

    def add_actor(self):
        self.lbf = self.make_actor()
        self.system.lb = self.lbf
        return self.lbf

    def clear_actors(self):
        self.system.lb = None

    @utx.skipIfMissingModules("espressomd.io.vtk")
    def test_vtk(self):
        """
        Check VTK files. Keep in mind the VTK module writes in single-precision.
        """
        dist = 1.5 * self.lattice.agrid
        actor = self.lbf
        actor.add_boundary_from_shape(
            espressomd.shapes.Wall(normal=[1, 0, 0], dist=dist))
        actor.add_boundary_from_shape(
            espressomd.shapes.Wall(normal=[-1, 0, 0], dist=-(self.system.box_l[0] - dist)))
        sphere = espressomd.shapes.Sphere(
            center=self.system.box_l // 2, radius=2.)
        actor.add_boundary_from_shape(sphere)

        n_steps = 4 if self.lb_params["gpu"] else 10
        shape = tuple(actor.shape)
        shape = (shape[0] - 4, *shape[1:])
        vtk_reader = espressomd.io.vtk.VTKReader()
        label_density = "density"
        label_velocity = "velocity_vector"
        label_pressure = "pressure_tensor"
        self.lbf[2, :, :].density = 1.3
        self.lbf[-3, :, :].density = 1.3

        with tempfile.TemporaryDirectory() as tmp_directory:
            root = pathlib.Path(tmp_directory)
            label_vtk_last_frame = f"test_vtk_{self.vtk_id}_last_frame"
            label_vtk_continuous = f"test_vtk_{self.vtk_id}_continuous"
            label_vtk_with_boundaries = f"test_vtk_{self.vtk_id}_with_boundaries"  # nopep8
            path_vtk_last_frame = root / label_vtk_last_frame / "simulation_step_0.vtu"
            path_vtk_continuous = [
                root / label_vtk_continuous / f"simulation_step_{i}.vtu" for i in range(n_steps)]
            path_vtk_with_boundaries = root / \
                label_vtk_with_boundaries / "simulation_step_0.vtu"
            filepaths = [path_vtk_last_frame] + path_vtk_continuous

            # write VTK files
            vtk_obs = self.write_obs()
            vtk_obj = self.vtk_class(
                identifier=label_vtk_continuous, delta_N=1, observables=vtk_obs,
                base_folder=root)
            actor.add_vtk_writer(vtk=vtk_obj)
            vtk_obj.disable()
            vtk_obj.enable()
            self.system.integrator.run(n_steps)
            vtk_obj = self.vtk_class(
                identifier=label_vtk_last_frame, delta_N=0, observables=vtk_obs,
                base_folder=root)
            actor.add_vtk_writer(vtk=vtk_obj)
            vtk_obj.write()
            if self.include_boundaries:
                # also write a snapshot that includes boundary cells and the
                # ``boundary`` mask observable
                vtk_obj_b = self.vtk_class(
                    identifier=label_vtk_with_boundaries, delta_N=0,
                    observables=self.write_obs(), base_folder=root,
                    include_boundaries=True, force_pvtu=True)
                actor.add_vtk_writer(vtk=vtk_obj_b)
                vtk_obj_b.write()
            self.assertEqual(sorted(vtk_obj.observables), sorted(vtk_obs))
            self.assertEqual(vtk_obj.valid_observables(), set(self.valid_obs))

            # check VTK files exist
            for filepath in filepaths:
                self.assertTrue(
                    filepath.exists(),
                    f"VTK file \"{filepath}\" not written to disk")
            for filepath in [path_vtk_last_frame.parent.with_suffix(".pvd"),
                             path_vtk_continuous[0].parent.with_suffix(".pvd")]:
                self.assertTrue(
                    filepath.exists(),
                    f"VTK summary file \"{filepath}\" not written to disk")

            # check velocity profile is symmetric at all time steps
            for filepath in filepaths:
                vtk_velocity = vtk_reader.parse(filepath)[label_velocity]
                v_profile = np.mean(
                    np.linalg.norm(vtk_velocity, axis=-1),
                    axis=(1, 2))
                np.testing.assert_allclose(
                    v_profile, v_profile[::-1], rtol=5e-5, atol=0.)

            # check scalar pressure is symmetric at all time steps
            for filepath in filepaths:
                vtk_pressure = vtk_reader.parse(filepath)[label_pressure]
                vtk_pressure = vtk_pressure.reshape(shape + (3, 3))
                p_profile = np.mean(
                    np.trace(vtk_pressure, axis1=-2, axis2=-1),
                    axis=(1, 2))
                np.testing.assert_allclose(
                    p_profile, p_profile[::-1], rtol=5e-5, atol=0.)

            # read VTK output of final time step
            last_frames = []
            for filepath in (path_vtk_last_frame, path_vtk_continuous[-1]):
                grids = vtk_reader.parse(filepath)
                last_frames.append((
                    grids[label_density],
                    grids[label_velocity],
                    grids[label_pressure].reshape(shape + (3, 3)),
                ))

            # check VTK output is identical in both continuous and manual mode
            for i in range(len(last_frames[0])):
                np.testing.assert_allclose(last_frames[0][i],
                                           last_frames[1][i], atol=1e-10)

            # build boundary mask and verify NaN only at boundaries
            inner_mask = np.copy(self.lbf[2:-2, :, :].is_boundary)
            for vtk_density, vtk_velocity, vtk_pressure in last_frames:
                nan_mask = np.isnan(vtk_density[:, :, :])
                np.testing.assert_array_equal(
                    nan_mask, inner_mask,
                    "NaN values in VTK data do not match boundary mask")
            lb_density = np.copy(self.lbf[2:-2, :, :].density)
            lb_velocity = np.copy(self.lbf[2:-2, :, :].velocity)
            lb_pressure = np.copy(self.lbf[2:-2, :, :].pressure_tensor)

            for vtk_density, vtk_velocity, vtk_pressure in last_frames:
                valid = ~np.isnan(vtk_density)
                np.testing.assert_allclose(
                    vtk_density[valid], lb_density[valid], rtol=1e-7, atol=0.)
                np.testing.assert_allclose(
                    vtk_velocity[valid], lb_velocity[valid], rtol=1e-7, atol=0.)
                np.testing.assert_allclose(
                    vtk_pressure[valid], lb_pressure[valid], rtol=1e-6, atol=0.)

            if self.include_boundaries:
                # check the include_boundaries snapshot: full lattice shape and
                # correct boundary mask in the two outer slabs
                full_shape = tuple(actor.shape)
                grids_b = vtk_reader.parse(path_vtk_with_boundaries)
                self.assertEqual(grids_b[label_density].shape, full_shape)
                self.assertEqual(grids_b["boundary"].shape, full_shape)
                expected_mask = np.zeros(full_shape, dtype=np.float32)
                expected_mask[:2, :, :] = 1.
                expected_mask[-2:, :, :] = 1.
                expected_mask[self.lattice.get_shape_bitmask(
                    shape=sphere)] = 1.
                np.testing.assert_array_equal(
                    grids_b["boundary"], expected_mask)
                np.testing.assert_array_equal(
                    np.asarray(actor[:, :, :].is_boundary, dtype=np.float32),
                    expected_mask)
                # the fluid region of the include_boundaries output matches the
                # filtered output
                np.testing.assert_allclose(
                    grids_b[label_density][2:-2, :, :], lb_density,
                    rtol=1e-7, atol=0.)
                # check that boundary cell values are written correctly
                # density and velocity in boundary region are available via the
                # boundary mask; the interior data matches the filtered output
                flat_boundary_mask = np.asarray(
                    grids_b["boundary"]).ravel().astype(bool)
                np.testing.assert_allclose(
                    grids_b[label_density].ravel()[flat_boundary_mask],
                    np.copy(self.actor[:, :, :].density).ravel()[
                        flat_boundary_mask],
                    rtol=1e-7, atol=0.)
                np.testing.assert_allclose(
                    grids_b[label_velocity].reshape(-1, 3)[flat_boundary_mask],
                    np.copy(
                        self.actor[:, :, :].velocity).reshape(-1, 3)[flat_boundary_mask],
                    rtol=1e-7, atol=0.)

    @utx.skipIfMissingModules("espressomd.io.vtk")
    def test_utf8_support(self):
        """Check UTF-8 support in filepaths and VTK identifiers."""
        with tempfile.TemporaryDirectory() as tmp_directory:
            root = pathlib.Path(tmp_directory) / "gemäß"
            label = "çåš"
            vtk_obs = self.write_obs()
            vtk_obj = self.vtk_class(
                identifier=label, delta_N=0, observables=vtk_obs, base_folder=root)
            self.lbf.add_vtk_writer(vtk=vtk_obj)
            self.assertEqual(vtk_obj.identifier, label)
            self.assertEqual(vtk_obj.base_folder, root)
            vtk_obj.write()
            path = root / label / "simulation_step_0.vtu"
            self.assertTrue(path.exists(), f"File \"{path}\" not found")


class TestEKVTK(TestVTK):
    include_boundaries = True

    valid_obs = ["density", "flux", "boundary"]
    valid_obs_poisson = ["potential"]

    def write_obs(self):
        obs = ["density", "flux"]
        if self.include_boundaries:
            obs.append("boundary")
        return obs

    def make_actor(self):
        return self.ek_class(
            lattice=self.lattice, density=1., diffusion=0.1, valency=0.1, kT=1.,
            advection=False, friction_coupling=False, tau=0.1, **self.ek_params)

    def add_actor(self):
        self.solver = self.ek_solver(
            lattice=self.lattice, permittivity=0.1, tau=0.1, **self.ek_params)
        self.species = self.make_actor()
        self.system.ekcontainer = espressomd.electrokinetics.EKContainer(
            tau=0.1, solver=self.solver)
        self.system.ekcontainer.add(self.species)
        return self.species

    def clear_actors(self):
        self.system.ekcontainer = None

    @utx.skipIfMissingModules("espressomd.io.vtk")
    def test_vtk(self):
        """
        Check VTK files. Keep in mind the VTK module writes in single-precision.
        """
        dist = 1.5 * self.lattice.agrid
        actor = self.species
        actor.add_boundary_from_shape(
            shape=espressomd.shapes.Wall(normal=[1, 0, 0], dist=dist),
            value=0.0, boundary_type=espressomd.electrokinetics.DensityBoundary)
        actor.add_boundary_from_shape(
            shape=espressomd.shapes.Wall(
                normal=[-1, 0, 0], dist=-(self.system.box_l[0] - dist)),
            value=0.0, boundary_type=espressomd.electrokinetics.DensityBoundary)
        sphere = espressomd.shapes.Sphere(
            center=self.system.box_l // 2, radius=2.)
        actor.add_boundary_from_shape(
            shape=sphere, value=0.0, boundary_type=espressomd.electrokinetics.DensityBoundary)
        actor[2, 0, 0].flux_boundary = espressomd.electrokinetics.FluxBoundary(
            [0.01, -0.01, 0.02])
        if isinstance(self.solver, espressomd.electrokinetics.EKNone):
            kx, ky, kz = np.pi / self.lattice.shape
            self.solver[:, :, :].potential = np.fromfunction(
                lambda i, j, k: np.cos(i * kx) *
                np.cos(j * ky) * np.cos(k * kz),
                self.lattice.shape, dtype=float)

        n_steps = 100
        shape = tuple(self.lattice.shape)
        shape = (shape[0] - 4, *shape[1:])
        vtk_reader = espressomd.io.vtk.VTKReader()
        label_density = "density"
        label_flux = "flux"
        label_potential = "potential"

        with tempfile.TemporaryDirectory() as tmp_directory:
            root = pathlib.Path(tmp_directory)
            label_vtk_last_frame = f"test_vtk_{self.vtk_id}_end"
            label_vtk_continuous = f"test_vtk_{self.vtk_id}_continuous"
            path_vtk_last_frame = root / label_vtk_last_frame / "simulation_step_0.vtu"
            path_vtk_continuous = [
                root / label_vtk_continuous / f"simulation_step_{i}.vtu" for i in range(n_steps)]
            filepaths = [path_vtk_last_frame] + path_vtk_continuous

            # write VTK files
            vtk_obs = self.write_obs()
            vtk_obj = self.vtk_class(
                identifier=label_vtk_continuous, delta_N=1, observables=vtk_obs,
                base_folder=root)
            actor.add_vtk_writer(vtk=vtk_obj)
            vtk_obj.disable()
            vtk_obj.enable()
            self.assertFalse(vtk_obj.force_pvtu)

            # prepare VTK Poisson
            label_vtk_poisson_last_frame = f"test_vtk_{self.vtk_id}_poisson_end"  # nopep8
            label_vtk_poisson_continuous = f"test_vtk_{self.vtk_id}_poisson_continuous"  # nopep8
            path_vtk_poisson_last_frame = root / \
                label_vtk_poisson_last_frame / "simulation_step_0.vti"
            path_vtk_poisson_continuous = [
                root / label_vtk_poisson_continuous / f"simulation_step_{i}.vtu" for i in range(n_steps)]
            filepaths_poisson = [
                path_vtk_poisson_last_frame] + path_vtk_poisson_continuous

            vtk_obs_poisson = list(self.valid_obs_poisson)
            vtk_obj_poisson = self.vtk_poisson_class(
                identifier=label_vtk_poisson_continuous, force_pvtu=True,
                observables=vtk_obs_poisson, base_folder=root, delta_N=1)
            self.assertTrue(vtk_obj_poisson.force_pvtu)
            self.solver.add_vtk_writer(vtk=vtk_obj_poisson)
            vtk_obj_poisson.disable()
            vtk_obj_poisson.enable()

            self.system.integrator.run(n_steps)

            # write manual files after integration
            vtk_obj = self.vtk_class(
                identifier=label_vtk_last_frame, delta_N=0, observables=vtk_obs,
                base_folder=root)
            actor.add_vtk_writer(vtk=vtk_obj)
            vtk_obj.write()
            self.assertEqual(sorted(vtk_obj.observables), sorted(vtk_obs))
            self.assertEqual(vtk_obj.valid_observables(), set(self.valid_obs))

            if self.include_boundaries:
                # also write a snapshot that includes boundary cells and the
                # ``boundary`` mask observable
                label_vtk_with_boundaries = f"test_vtk_{self.vtk_id}_with_boundaries"  # nopep8
                path_vtk_with_boundaries = root / \
                    label_vtk_with_boundaries / "simulation_step_0.vtu"
                vtk_obj_b = self.vtk_class(
                    identifier=label_vtk_with_boundaries, delta_N=0,
                    observables=self.write_obs(), base_folder=root,
                    include_boundaries=True, force_pvtu=True)
                actor.add_vtk_writer(vtk=vtk_obj_b)
                vtk_obj_b.write()

            vtk_obj_poisson = self.vtk_poisson_class(
                identifier=label_vtk_poisson_last_frame, force_pvtu=False,
                observables=vtk_obs_poisson, base_folder=root, delta_N=0)
            self.assertFalse(vtk_obj_poisson.force_pvtu)
            self.solver.add_vtk_writer(vtk=vtk_obj_poisson)
            vtk_obj_poisson.write()
            self.assertEqual(sorted(vtk_obj_poisson.observables),
                             sorted(vtk_obs_poisson))
            self.assertEqual(vtk_obj_poisson.valid_observables(),
                             set(self.valid_obs_poisson))

            # check VTK files exist
            for filepath in filepaths + filepaths_poisson:
                self.assertTrue(
                    filepath.exists(),
                    f"VTK file \"{filepath}\" not written to disk")
            for filepath in [path_vtk_last_frame.parent.with_suffix(".pvd"),
                             path_vtk_continuous[0].parent.with_suffix(".pvd"),
                             path_vtk_poisson_last_frame.parent.with_suffix(
                                 ".pvd"),
                             path_vtk_poisson_continuous[0].parent.with_suffix(".pvd")]:
                self.assertTrue(
                    filepath.exists(),
                    f"VTK summary file \"{filepath}\" not written to disk")

            # read VTK output of final time step
            last_frames = []
            last_frames_flux = []
            for filepath in (path_vtk_last_frame, path_vtk_continuous[-1],):
                grids = vtk_reader.parse(filepath)
                last_frames.append(grids[label_density])
                last_frames_flux.append(grids[label_flux])

            last_frames_poisson = []
            for filepath in (path_vtk_poisson_last_frame,
                             path_vtk_poisson_continuous[-1],):
                grids = vtk_reader.parse(filepath)
                last_frames_poisson.append(grids[label_potential])

            # check VTK output is identical in both continuous and manual mode
            for i in range(len(last_frames[0])):
                np.testing.assert_allclose(last_frames[0][i],
                                           last_frames[1][i], atol=1e-10)
                np.testing.assert_allclose(last_frames_flux[0][i],
                                           last_frames_flux[1][i], atol=1e-10)

            for i in range(len(last_frames_poisson[0])):
                np.testing.assert_allclose(last_frames_poisson[0][i],
                                           last_frames_poisson[1][i], atol=1e-10)

            # check VTK values match node values in the final time step
            tol = {"rtol": 5e-7, "atol": 1e-12}
            ek_inner_mask = np.copy(actor[2:-2, :, :].is_boundary)

            for vtk_density in last_frames:
                nan_mask = np.isnan(vtk_density[:, :, :])
                np.testing.assert_array_equal(
                    nan_mask, ek_inner_mask,
                    "NaN values in VTK density do not match boundary mask")

            ek_density = np.copy(actor[2:-2, :, :].density)

            for vtk_density in last_frames:
                valid = ~np.isnan(vtk_density)
                np.testing.assert_allclose(
                    vtk_density[valid], ek_density[valid], **tol)

            ek_flux = np.copy(actor[2:-2, :, :].flux)
            ek_flux_mask = np.repeat(
                ek_inner_mask[..., np.newaxis], 3, axis=-1)
            for vtk_flux in last_frames_flux:
                nan_mask = np.isnan(vtk_flux[:, :, :])
                np.testing.assert_array_equal(
                    nan_mask, ek_flux_mask,
                    "NaN values in VTK flux do not match boundary mask")
            for vtk_flux in last_frames_flux:
                valid = ~np.isnan(vtk_flux)
                np.testing.assert_allclose(
                    vtk_flux[valid], ek_flux[valid], **tol)

            ek_potential = np.copy(self.solver[:, :, :].potential)

            for vtk_potential in last_frames_poisson:
                valid = ~np.isnan(vtk_potential)
                np.testing.assert_allclose(
                    vtk_potential[valid], ek_potential[valid], **tol)

            if self.include_boundaries:
                # check the include_boundaries snapshot: full lattice shape and
                # correct boundary mask in the two outer slabs
                full_shape = tuple(self.lattice.shape)
                grids_b = vtk_reader.parse(path_vtk_with_boundaries)
                self.assertEqual(grids_b[label_density].shape, full_shape)
                self.assertEqual(grids_b["boundary"].shape, full_shape)
                expected_mask = np.zeros(full_shape, dtype=np.float32)
                expected_mask[:2, :, :] = 1.
                expected_mask[-2:, :, :] = 1.
                expected_mask[self.lattice.get_shape_bitmask(
                    shape=sphere)] = 1.
                np.testing.assert_array_equal(
                    grids_b["boundary"], expected_mask)
                np.testing.assert_allclose(
                    grids_b[label_density][2:-2, :, :], ek_density, **tol)
                # check that boundary cell values are written correctly
                boundary_mask = actor[:, :, :].is_boundary
                vtk_boundary_density = grids_b[label_density]
                vtk_boundary_flux = grids_b[label_flux]
                np.testing.assert_allclose(
                    vtk_boundary_density[boundary_mask],
                    np.copy(self.species[:, :, :].density)[boundary_mask],
                    **tol)
                np.testing.assert_allclose(
                    vtk_boundary_flux[boundary_mask],
                    np.copy(self.species[:, :, :].flux)[boundary_mask],
                    **tol)

        expected_writers = 3 if self.include_boundaries else 2
        self.assertEqual(len(actor.vtk_writers), expected_writers)
        actor.clear_vtk_writers()
        self.assertEqual(len(actor.vtk_writers), 0)


@utx.skipIfMissingFeatures(["WALBERLA"])
class LBWalberlaVTKDoublePrecisionCPU(TestLBVTK, ut.TestCase):
    vtk_class = espressomd.lb.VTKOutput
    lattice_class = espressomd.lb.Lattice
    lb_class = espressomd.lb.LBFluid
    lb_params = {"single_precision": False, "gpu": False}
    vtk_id = "lb_double_precision_cpu"


@utx.skipIfMissingFeatures(["WALBERLA"])
class LBWalberlaVTKDoublePrecisionBlocksCPU(TestLBVTK, ut.TestCase):
    vtk_class = espressomd.lb.VTKOutput
    lattice_class = espressomd.lb.Lattice
    lb_class = espressomd.lb.LBFluid
    lb_params = {"single_precision": False, "gpu": False}
    # use more than one waLBerla block per MPI rank: with box_l=[6,7,3] and
    # agrid=0.5 the lattice is 12x14x6, so [2,1,1] splits it into two blocks
    # along x on a single rank. The VTK field writers must export each block's
    # own data, not the last block's data for every block.
    lattice_params = {"blocks_per_mpi_rank": [2, 1, 1]}
    vtk_id = "lb_double_precision_blocks_cpu"


@utx.skipIfMissingGPU()
@utx.skipIfMissingFeatures(["WALBERLA", "CUDA"])
class LBWalberlaVTKSinglePrecisionGPU(TestLBVTK, ut.TestCase):
    vtk_class = espressomd.lb.VTKOutput
    lattice_class = espressomd.lb.Lattice
    lb_class = espressomd.lb.LBFluid
    lb_params = {"single_precision": True, "gpu": True}
    vtk_id = "lb_single_precision_gpu"


@utx.skipIfMissingFeatures(["WALBERLA", "WALBERLA_FFT"])
class EKWalberlaVTKDoublePrecisionCPU(TestEKVTK, ut.TestCase):
    vtk_class = espressomd.electrokinetics.VTKOutput
    vtk_poisson_class = espressomd.electrokinetics.VTKPoissonOutput
    lattice_class = espressomd.electrokinetics.Lattice
    ek_class = espressomd.electrokinetics.EKSpecies
    ek_solver = espressomd.electrokinetics.EKFFT
    ek_params = {"single_precision": False, "gpu": False}
    vtk_id = "ek_double_precision_cpu"


@utx.skipIfMissingGPU()
@utx.skipIfMissingFeatures(["WALBERLA", "WALBERLA_FFT", "CUDA"])
class EKWalberlaVTKSinglePrecisionGPU(TestEKVTK, ut.TestCase):
    vtk_class = espressomd.electrokinetics.VTKOutput
    vtk_poisson_class = espressomd.electrokinetics.VTKPoissonOutput
    lattice_class = espressomd.electrokinetics.Lattice
    ek_class = espressomd.electrokinetics.EKSpecies
    ek_solver = espressomd.electrokinetics.EKFFT
    ek_params = {"single_precision": True, "gpu": True}
    vtk_id = "ek_single_precision_gpu"


@utx.skipIfMissingFeatures(["WALBERLA"])
class EKWalberlaVTKDoublePrecisionEKNoneCPU(TestEKVTK, ut.TestCase):
    vtk_class = espressomd.electrokinetics.VTKOutput
    vtk_poisson_class = espressomd.electrokinetics.VTKPoissonOutput
    lattice_class = espressomd.electrokinetics.Lattice
    ek_class = espressomd.electrokinetics.EKSpecies
    ek_solver = espressomd.electrokinetics.EKNone
    ek_params = {"single_precision": False, "gpu": False}
    vtk_id = "ek_double_precision_cpu_eknone"


@utx.skipIfMissingGPU()
@utx.skipIfMissingFeatures(["WALBERLA", "CUDA"])
class EKWalberlaVTKSinglePrecisionEKNoneGPU(TestEKVTK, ut.TestCase):
    vtk_class = espressomd.electrokinetics.VTKOutput
    vtk_poisson_class = espressomd.electrokinetics.VTKPoissonOutput
    lattice_class = espressomd.electrokinetics.Lattice
    ek_class = espressomd.electrokinetics.EKSpecies
    ek_solver = espressomd.electrokinetics.EKNone
    ek_params = {"single_precision": True, "gpu": True}
    vtk_id = "ek_single_precision_gpu_eknone"


@utx.skipIfMissingGPU()
@utx.skipIfMissingFeatures(["WALBERLA", "CUDA"])
class LBWalberlaVTKDoublePrecisionGPU_NoBoundaries(TestLBVTK, ut.TestCase):
    vtk_class = espressomd.lb.VTKOutput
    lattice_class = espressomd.lb.Lattice
    lb_class = espressomd.lb.LBFluid
    lb_params = {"single_precision": False, "gpu": True}
    vtk_id = "lb_double_precision_gpu_no_boundaries"
    include_boundaries = False


@utx.skipIfMissingFeatures(["WALBERLA"])
class LBWalberlaVTKSinglePrecisionCPU_NoBoundaries(TestLBVTK, ut.TestCase):
    vtk_class = espressomd.lb.VTKOutput
    lattice_class = espressomd.lb.Lattice
    lb_class = espressomd.lb.LBFluid
    lb_params = {"single_precision": True, "gpu": False}
    vtk_id = "lb_single_precision_cpu_no_boundaries"
    include_boundaries = False


@utx.skipIfMissingGPU()
@utx.skipIfMissingFeatures(["WALBERLA", "WALBERLA_FFT", "CUDA"])
class EKWalberlaVTKDoublePrecisionGPU_NoBoundaries(TestEKVTK, ut.TestCase):
    vtk_class = espressomd.electrokinetics.VTKOutput
    vtk_poisson_class = espressomd.electrokinetics.VTKPoissonOutput
    lattice_class = espressomd.electrokinetics.Lattice
    ek_class = espressomd.electrokinetics.EKSpecies
    ek_solver = espressomd.electrokinetics.EKFFT
    ek_params = {"single_precision": False, "gpu": True}
    vtk_id = "ek_double_precision_gpu_no_boundaries"
    include_boundaries = False


@utx.skipIfMissingFeatures(["WALBERLA", "WALBERLA_FFT"])
class EKWalberlaVTKSinglePrecisionCPU_NoBoundaries(TestEKVTK, ut.TestCase):
    vtk_class = espressomd.electrokinetics.VTKOutput
    vtk_poisson_class = espressomd.electrokinetics.VTKPoissonOutput
    lattice_class = espressomd.electrokinetics.Lattice
    ek_class = espressomd.electrokinetics.EKSpecies
    ek_solver = espressomd.electrokinetics.EKFFT
    ek_params = {"single_precision": True, "gpu": False}
    vtk_id = "ek_single_precision_cpu_no_boundaries"
    include_boundaries = False


if __name__ == "__main__":
    ut.main()
