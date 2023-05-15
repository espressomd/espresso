#
# Copyright (C) 2010-2023 The ESPResSo project
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
    system = espressomd.System(box_l=[6, 7, 5])
    system.time_step = 0.1
    system.cell_system.skin = 0.4

    def setUp(self):
        self.lattice = self.lattice_class(n_ghost_layers=1, agrid=0.5)
        self.actor = self.add_actor()

    def tearDown(self):
        self.clear_actors()

    @ut.skipIf(system.cell_system.get_state()["n_nodes"] > 4,
               "this test is slow on more than 4 MPI ranks")
    def test_exceptions(self):
        label_invalid_obs = f"test_vtk_{self.vtk_id}_invalid_obs"
        error_msg = rf"Only the following VTK observables are supported: \[{repr(sorted(self.valid_obs))[1:-1]}\], got 'dens'"
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


class TestLBVTK(TestVTK):

    valid_obs = ["density", "velocity_vector", "pressure_tensor"]

    def make_actor(self):
        return self.lb_class(
            lattice=self.lattice, tau=0.1, density=1., kinematic_viscosity=1.,
            ext_force_density=[0., 0.03, 0.], **self.lb_params)

    def add_actor(self):
        self.lbf = self.make_actor()
        self.system.actors.add(self.lbf)
        return self.lbf

    def clear_actors(self):
        self.system.actors.clear()

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

        n_steps = 10
        shape = tuple(actor.shape)
        shape = (shape[0] - 4, *shape[1:])
        vtk_reader = espressomd.io.vtk.VTKReader()
        label_density = "density"
        label_velocity = "velocity_vector"
        label_pressure = "pressure_tensor"

        with tempfile.TemporaryDirectory() as tmp_directory:
            path_vtk_root = pathlib.Path(tmp_directory)
            label_vtk_end = f"test_vtk_{self.vtk_id}_end"
            label_vtk_continuous = f"test_vtk_{self.vtk_id}_continuous"
            path_vtk_end = path_vtk_root / label_vtk_end / "simulation_step_0.vtu"
            path_vtk_continuous = [
                path_vtk_root / label_vtk_continuous / f"simulation_step_{i}.vtu" for i in range(n_steps)]
            filepaths = [path_vtk_end] + path_vtk_continuous

            # write VTK files
            vtk_obs = list(self.valid_obs)
            vtk_obj = self.vtk_class(
                identifier=label_vtk_continuous, delta_N=1, observables=vtk_obs,
                base_folder=str(path_vtk_root))
            actor.add_vtk_writer(vtk=vtk_obj)
            vtk_obj.disable()
            vtk_obj.enable()
            self.system.integrator.run(n_steps)
            vtk_obj = self.vtk_class(
                identifier=label_vtk_end, delta_N=0, observables=vtk_obs,
                base_folder=str(path_vtk_root))
            actor.add_vtk_writer(vtk=vtk_obj)
            vtk_obj.write()
            self.assertEqual(sorted(vtk_obj.observables), sorted(vtk_obs))
            self.assertEqual(vtk_obj.valid_observables(), set(self.valid_obs))

            # check VTK files exist
            for filepath in filepaths:
                self.assertTrue(
                    filepath.exists(),
                    f"VTK file \"{filepath}\" not written to disk")
            for filepath in [path_vtk_end.parent.with_suffix(".pvd"),
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
            for filepath in (path_vtk_end, path_vtk_continuous[-1]):
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

            # check VTK values match node values in the final time step
            lb_density = np.copy(self.lbf[2:-2, :, :].density)
            lb_velocity = np.copy(self.lbf[2:-2, :, :].velocity)
            lb_pressure = np.copy(self.lbf[2:-2, :, :].pressure_tensor)

            for vtk_density, vtk_velocity, vtk_pressure in last_frames:
                np.testing.assert_allclose(
                    vtk_density, lb_density, rtol=1e-10, atol=0.)
                np.testing.assert_allclose(
                    vtk_velocity, lb_velocity, rtol=1e-7, atol=0.)
                np.testing.assert_allclose(
                    vtk_pressure, lb_pressure, rtol=1e-6, atol=0.)


class TestEKVTK(TestVTK):

    valid_obs = ["density"]

    def make_actor(self):
        return self.ek_class(
            lattice=self.lattice, density=1., diffusion=0.1, valency=0.,
            advection=False, friction_coupling=False, tau=0.1, **self.ek_params)

    def add_actor(self):
        self.solver = self.ek_solver(lattice=self.lattice)
        self.species = self.make_actor()
        self.system.ekcontainer.tau = 0.1
        self.system.ekcontainer.solver = self.solver
        self.system.ekcontainer.add(self.species)
        return self.species

    def clear_actors(self):
        self.system.ekcontainer.clear()

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

        n_steps = 100
        shape = tuple(self.lattice.shape)
        shape = (shape[0] - 4, *shape[1:])
        vtk_reader = espressomd.io.vtk.VTKReader()
        label_density = "density"

        with tempfile.TemporaryDirectory() as tmp_directory:
            path_vtk_root = pathlib.Path(tmp_directory)
            label_vtk_end = f"test_vtk_{self.vtk_id}_end"
            label_vtk_continuous = f"test_vtk_{self.vtk_id}_continuous"
            path_vtk_end = path_vtk_root / label_vtk_end / "simulation_step_0.vtu"
            path_vtk_continuous = [
                path_vtk_root / label_vtk_continuous / f"simulation_step_{i}.vtu" for i in range(n_steps)]
            filepaths = [path_vtk_end] + path_vtk_continuous

            # write VTK files
            vtk_obs = list(self.valid_obs)
            vtk_obj = self.vtk_class(
                identifier=label_vtk_continuous, delta_N=1, observables=vtk_obs,
                base_folder=str(path_vtk_root))
            actor.add_vtk_writer(vtk=vtk_obj)
            vtk_obj.disable()
            vtk_obj.enable()
            self.system.integrator.run(n_steps)
            vtk_obj = self.vtk_class(
                identifier=label_vtk_end, delta_N=0, observables=vtk_obs,
                base_folder=str(path_vtk_root))
            actor.add_vtk_writer(vtk=vtk_obj)
            vtk_obj.write()
            self.assertEqual(sorted(vtk_obj.observables), sorted(vtk_obs))
            self.assertEqual(vtk_obj.valid_observables(), set(self.valid_obs))

            # check VTK files exist
            for filepath in filepaths:
                self.assertTrue(
                    filepath.exists(),
                    f"VTK file \"{filepath}\" not written to disk")
            for filepath in [path_vtk_end.parent.with_suffix(".pvd"),
                             path_vtk_continuous[0].parent.with_suffix(".pvd")]:
                self.assertTrue(
                    filepath.exists(),
                    f"VTK summary file \"{filepath}\" not written to disk")

            # read VTK output of final time step
            last_frames = []
            for filepath in (path_vtk_end, path_vtk_continuous[-1]):
                grids = vtk_reader.parse(filepath)
                last_frames.append(grids[label_density])

            # check VTK output is identical in both continuous and manual mode
            for i in range(len(last_frames[0])):
                np.testing.assert_allclose(last_frames[0][i],
                                           last_frames[1][i], atol=1e-10)

            # check VTK values match node values in the final time step
            ek_density = np.copy(actor[2:-2, :, :].density)

            for vtk_density in last_frames:
                np.testing.assert_allclose(
                    vtk_density, ek_density, rtol=5e-7)

        self.assertEqual(len(actor.vtk_writers), 2)
        actor.clear_vtk_writers()
        self.assertEqual(len(actor.vtk_writers), 0)


@utx.skipIfMissingFeatures("WALBERLA")
class LBWalberlaWrite(TestLBVTK, ut.TestCase):
    vtk_class = espressomd.lb.VTKOutput
    lattice_class = espressomd.lb.LatticeWalberla
    lb_class = espressomd.lb.LBFluidWalberla
    lb_params = {"single_precision": False}
    vtk_id = "lb_double_precision"


@utx.skipIfMissingFeatures("WALBERLA")
class LBWalberlaWriteSinglePrecision(TestLBVTK, ut.TestCase):
    vtk_class = espressomd.lb.VTKOutput
    lattice_class = espressomd.lb.LatticeWalberla
    lb_class = espressomd.lb.LBFluidWalberla
    lb_params = {"single_precision": True}
    vtk_id = "lb_single_precision"


@utx.skipIfMissingFeatures("WALBERLA")
class EKWalberlaWrite(TestEKVTK, ut.TestCase):
    vtk_class = espressomd.electrokinetics.VTKOutput
    lattice_class = espressomd.electrokinetics.LatticeWalberla
    ek_class = espressomd.electrokinetics.EKSpecies
    ek_solver = espressomd.electrokinetics.EKNone
    ek_params = {"single_precision": False}
    vtk_id = "ek_double_precision"


@utx.skipIfMissingFeatures("WALBERLA")
class EKWalberlaWriteSinglePrecision(TestEKVTK, ut.TestCase):
    vtk_class = espressomd.electrokinetics.VTKOutput
    lattice_class = espressomd.electrokinetics.LatticeWalberla
    ek_class = espressomd.electrokinetics.EKSpecies
    ek_solver = espressomd.electrokinetics.EKNone
    ek_params = {"single_precision": True}
    vtk_id = "ek_single_precision"


if __name__ == "__main__":
    ut.main()
