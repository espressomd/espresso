#
# Copyright (C) 2010-2022 The ESPResSo project
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

# pylint: disable=undefined-variable
import unittest as ut
import unittest_decorators as utx
import unittest_generator as utg
import numpy as np
import contextlib
import pathlib
import sys

import espressomd
import espressomd.checkpointing
import espressomd.electrostatics
import espressomd.magnetostatics
import espressomd.io.writer  # pylint: disable=unused-import
import espressomd.lees_edwards
import espressomd.virtual_sites
import espressomd.integrate
import espressomd.shapes
import espressomd.constraints
import espressomd.lb
import espressomd.electrokinetics
with contextlib.suppress(ImportError):
    import espressomd.io.vtk

with contextlib.suppress(ImportError):
    import h5py  # h5py has to be imported *after* espressomd (MPI)

config = utg.TestGenerator()
is_gpu_available = espressomd.gpu_available()
modes = config.get_modes()
has_lb_mode = ('LB.WALBERLA' in modes and espressomd.has_features('WALBERLA')
               and ('LB.CPU' in modes or 'LB.GPU' in modes and is_gpu_available))
has_p3m_mode = 'P3M.CPU' in modes or 'P3M.GPU' in modes and is_gpu_available


class CheckpointTest(ut.TestCase):

    checkpoint = espressomd.checkpointing.Checkpoint(
        **config.get_checkpoint_params())
    checkpoint.load(0)
    path_cpt_root = pathlib.Path(checkpoint.checkpoint_dir)

    @classmethod
    def setUpClass(cls):
        cls.ref_box_l = np.array([12.0, 8.0, 16.0])
        if 'DP3M' in modes:
            cls.ref_box_l = np.array([16.0, 16.0, 16.0])
        cls.ref_periodicity = np.array([True, True, True])
        if espressomd.has_features('STOKESIAN_DYNAMICS') and (
                'THERM.SDM' in modes or 'INT.SDM' in modes):
            cls.ref_periodicity = np.array([False, False, False])

    def get_active_actor_of_type(self, actor_type):
        for actor in system.actors.active_actors:
            if isinstance(actor, actor_type):
                return actor
        self.fail(
            f"system doesn't have an actor of type {actor_type.__name__}")

    def test_get_active_actor_of_type(self):
        if system.actors.active_actors:
            actor = system.actors.active_actors[0]
            self.assertEqual(self.get_active_actor_of_type(type(actor)), actor)
        with self.assertRaisesRegex(AssertionError, "system doesn't have an actor of type Wall"):
            self.get_active_actor_of_type(espressomd.shapes.Wall)

    @utx.skipIfMissingFeatures(["WALBERLA"])
    @ut.skipIf(not has_lb_mode, "Skipping test due to missing LB feature.")
    def test_lb_fluid(self):
        lbf = self.get_active_actor_of_type(espressomd.lb.LBFluidWalberla)
        cpt_mode = 0 if 'LB.ASCII' in modes else 1
        cpt_root = pathlib.Path(self.checkpoint.checkpoint_dir)
        cpt_path = str(cpt_root / "lb") + "{}.cpt"

        # LB boundaries are loaded at the same time as LB populations
        np.testing.assert_equal(np.copy(lbf[:, :, :].velocity), 0.)
        np.testing.assert_equal(
            np.copy(lbf[:, :, :].is_boundary.astype(int)), 0)

        # check exception mechanism with corrupted LB checkpoint files
        with self.assertRaisesRegex(RuntimeError, 'EOF found'):
            lbf.load_checkpoint(cpt_path.format("-missing-data"), cpt_mode)
        with self.assertRaisesRegex(RuntimeError, 'extra data found, expected EOF'):
            lbf.load_checkpoint(cpt_path.format("-extra-data"), cpt_mode)
        if cpt_mode == 0:
            with self.assertRaisesRegex(RuntimeError, 'incorrectly formatted data'):
                lbf.load_checkpoint(cpt_path.format("-wrong-format"), cpt_mode)
            with self.assertRaisesRegex(RuntimeError, 'grid dimensions mismatch'):
                lbf.load_checkpoint(cpt_path.format("-wrong-boxdim"), cpt_mode)
            with self.assertRaisesRegex(RuntimeError, 'population size mismatch'):
                lbf.load_checkpoint(
                    cpt_path.format("-wrong-popsize"), cpt_mode)
        with self.assertRaisesRegex(RuntimeError, 'could not open file'):
            lbf.load_checkpoint(cpt_path.format("-unknown"), cpt_mode)

        # load the valid LB checkpoint file
        lbf.load_checkpoint(cpt_path.format(""), cpt_mode)
        precision = 8 if "LB.WALBERLA" in modes else 5
        m = np.pi / 12
        nx = lbf.shape[0]
        ny = lbf.shape[1]
        nz = lbf.shape[2]
        grid_3D = np.fromfunction(
            lambda i, j, k: np.cos(i * m) * np.cos(j * m) * np.cos(k * m),
            (nx, ny, nz), dtype=float)
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    np.testing.assert_almost_equal(
                        np.copy(lbf[i, j, k].population),
                        grid_3D[i, j, k] * np.arange(1, 20),
                        decimal=precision)
                    np.testing.assert_almost_equal(
                        np.copy(lbf[i, j, k].last_applied_force),
                        grid_3D[i, j, k] * np.arange(1, 4),
                        decimal=precision)
        state = lbf.get_params()
        reference = {
            "agrid": 2.0,
            "kinematic_viscosity": 1.3,
            "density": 1.5,
            "tau": 0.01}
        for key in reference:
            self.assertIn(key, state)
            np.testing.assert_allclose(np.copy(state[key]), reference[key],
                                       atol=1E-7, err_msg=f"{key} differs")
        self.assertTrue(lbf.is_active)
        self.assertFalse(lbf.single_precision)

        # check boundary objects
        slip_velocity1 = np.array([1e-4, 1e-4, 0.])
        slip_velocity2 = np.array([0., 0., 0.])
        # check boundary flag
        np.testing.assert_equal(
            np.copy(lbf[0, :, :].is_boundary.astype(int)), 1)
        np.testing.assert_equal(
            np.copy(lbf[-1, :, :].is_boundary.astype(int)), 1)
        np.testing.assert_equal(
            np.copy(lbf[1:-1, :, :].is_boundary.astype(int)), 0)
        # check boundary conditions
        for node in lbf[0, :, :]:
            np.testing.assert_allclose(np.copy(node.velocity), slip_velocity1)
        for node in lbf[-1, :, :]:
            np.testing.assert_allclose(np.copy(node.velocity), slip_velocity2)
        for node in lbf[2, :, :]:
            np.testing.assert_allclose(np.copy(node.velocity), 0.)
        # remove boundaries
        lbf.clear_boundaries()
        np.testing.assert_equal(
            np.copy(lbf[:, :, :].is_boundary.astype(int)), 0)

    @utx.skipIfMissingFeatures(["WALBERLA"])
    @ut.skipIf(not has_lb_mode, "Skipping test due to missing EK feature.")
    def test_ek_species(self):
        cpt_mode = 0 if 'LB.ASCII' in modes else 1
        cpt_root = pathlib.Path(self.checkpoint.checkpoint_dir)
        cpt_path = str(cpt_root / "ek") + "{}.cpt"

        self.assertEqual(len(system.ekcontainer), 1)
        ek_species = system.ekcontainer[0]
        self.assertTrue(
            system.ekcontainer.call_method("is_poisson_solver_set"))
        self.assertAlmostEqual(system.ekcontainer.tau, system.time_step,
                               delta=1e-7)
        self.assertIsInstance(system.ekcontainer.solver,
                              espressomd.electrokinetics.EKNone)

        # check exception mechanism with corrupted LB checkpoint files
        with self.assertRaisesRegex(RuntimeError, 'EOF found'):
            ek_species.load_checkpoint(
                cpt_path.format("-missing-data"), cpt_mode)
        with self.assertRaisesRegex(RuntimeError, 'extra data found, expected EOF'):
            ek_species.load_checkpoint(
                cpt_path.format("-extra-data"), cpt_mode)
        if cpt_mode == 0:
            with self.assertRaisesRegex(RuntimeError, 'incorrectly formatted data'):
                ek_species.load_checkpoint(
                    cpt_path.format("-wrong-format"), cpt_mode)
            with self.assertRaisesRegex(RuntimeError, 'grid dimensions mismatch'):
                ek_species.load_checkpoint(
                    cpt_path.format("-wrong-boxdim"), cpt_mode)
        with self.assertRaisesRegex(RuntimeError, 'could not open file'):
            ek_species.load_checkpoint(cpt_path.format("-unknown"), cpt_mode)

        ek_species.load_checkpoint(cpt_path.format(""), cpt_mode)

        precision = 8 if "LB.WALBERLA" in modes else 5
        m = np.pi / 12
        nx = ek_species.lattice.shape[0]
        ny = ek_species.lattice.shape[1]
        nz = ek_species.lattice.shape[2]
        grid_3D = np.fromfunction(
            lambda i, j, k: np.cos(i * m) * np.cos(j * m) * np.cos(k * m),
            (nx, ny, nz), dtype=float)
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    np.testing.assert_almost_equal(
                        np.copy(ek_species[i, j, k].density),
                        grid_3D[i, j, k], decimal=precision)

        state = ek_species.get_params()
        reference = {
            "density": 1.5,
            "diffusion": 0.2,
            "kT": 2.0,
            "valency": 0.1,
            "ext_efield": [0.1, 0.2, 0.3],
            "advection": False,
            "friction_coupling": False,
            "tau": 0.01}
        for key in reference:
            self.assertIn(key, state)
            np.testing.assert_allclose(np.copy(state[key]), reference[key],
                                       atol=1E-7, err_msg=f"{key} differs")
        # self.assertFalse(ek_species.is_active)
        self.assertFalse(ek_species.single_precision)

        def generator(value, shape):
            value_grid = np.tile(value, shape)
            if value_grid.shape[-1] == 1:
                value_grid = np.squeeze(value_grid, axis=-1)
            return value_grid

        # check boundary objects
        dens1 = 1.
        dens2 = 2.
        flux1 = 1e-3 * np.array([1., 2., 3.])
        flux2 = 1e-3 * np.array([4., 5., 6.])
        boundaries = [("density", dens1, dens2), ("flux", flux1, flux2)]
        for attr, value1, value2 in boundaries:
            accessor = np.vectorize(
                lambda obj: np.copy(getattr(obj, attr)),
                signature=f"()->({'n' if attr == 'flux' else ''})")
            slice1 = ek_species[0, :, :]
            slice2 = ek_species[-1, :, :]
            slice3 = ek_species[1:-1, :, :]
            # check boundary flag

            np.testing.assert_equal(np.copy(slice1.is_boundary), True)
            np.testing.assert_equal(np.copy(slice2.is_boundary), True)
            np.testing.assert_equal(np.copy(slice3.is_boundary), False)
            # check boundary conditions
            field = f"{attr}_boundary"
            shape = list(ek_species.shape)[-2:] + [1]
            np.testing.assert_allclose(
                accessor(np.copy(getattr(slice1, field))),
                generator(value1, shape))
            np.testing.assert_allclose(
                accessor(np.copy(getattr(slice2, field))),
                generator(value2, shape))

        ek_species.clear_density_boundaries()
        ek_species.clear_flux_boundaries()
        np.testing.assert_equal(
            np.copy(ek_species[:, :, :].is_boundary), False)

    @utx.skipIfMissingFeatures(["WALBERLA"])
    @ut.skipIf(not has_lb_mode, "Skipping test due to missing LB feature.")
    def test_lb_vtk(self):
        lbf = self.get_active_actor_of_type(espressomd.lb.LBFluidWalberla)
        self.assertEqual(len(lbf.vtk_writers), 2)
        vtk_suffix = config.test_name
        key_auto = f"vtk_out/auto_lb_{vtk_suffix}"
        vtk_auto = lbf.vtk_writers[0]
        self.assertIsInstance(vtk_auto, espressomd.lb.VTKOutput)
        self.assertEqual(vtk_auto.vtk_uid, key_auto)
        self.assertEqual(vtk_auto.delta_N, 1)
        self.assertFalse(vtk_auto.enabled)
        self.assertEqual(set(vtk_auto.observables),
                         {"density", "velocity_vector"})
        self.assertIn(
            f"write to '{key_auto}' every 1 LB steps (disabled)>", repr(vtk_auto))
        key_manual = f"vtk_out/manual_lb_{vtk_suffix}"
        vtk_manual = lbf.vtk_writers[1]
        self.assertIsInstance(vtk_manual, espressomd.lb.VTKOutput)
        self.assertEqual(vtk_manual.vtk_uid, key_manual)
        self.assertEqual(vtk_manual.delta_N, 0)
        self.assertEqual(set(vtk_manual.observables), {"density"})
        self.assertIn(f"write to '{key_manual}' on demand>", repr(vtk_manual))
        # check file numbering when resuming VTK write operations
        vtk_root = pathlib.Path("vtk_out") / f"manual_lb_{vtk_suffix}"
        filename = "simulation_step_{}.vtu"
        self.assertTrue((vtk_root / filename.format(0)).exists())
        self.assertFalse((vtk_root / filename.format(1)).exists())
        self.assertFalse((vtk_root / filename.format(2)).exists())
        # check VTK objects are still synchronized with their LB objects
        old_density = lbf[0, 0, 0].density
        new_density = 1.5 * old_density
        lbf[0, 0, 0].density = new_density
        vtk_manual.write()
        lbf[0, 0, 0].density = old_density
        self.assertTrue((vtk_root / filename.format(0)).exists())
        self.assertTrue((vtk_root / filename.format(1)).exists())
        self.assertFalse((vtk_root / filename.format(2)).exists())
        if "espressomd.io.vtk" in sys.modules:
            vtk_reader = espressomd.io.vtk.VTKReader()
            vtk_data = vtk_reader.parse(vtk_root / filename.format(1))
            lb_density = vtk_data["density"]
            self.assertAlmostEqual(
                lb_density[0, 0, 0], new_density, delta=1e-5)

    @utx.skipIfMissingFeatures(["WALBERLA"])
    @ut.skipIf(not has_lb_mode, "Skipping test due to missing EK feature.")
    def test_ek_vtk(self):
        ek_species = system.ekcontainer[0]
        vtk_suffix = config.test_name
        key_auto = f"vtk_out/auto_ek_{vtk_suffix}"
        vtk_auto = ek_species.vtk_writers[0]
        self.assertIsInstance(vtk_auto, espressomd.electrokinetics.VTKOutput)
        self.assertEqual(vtk_auto.vtk_uid, key_auto)
        self.assertEqual(vtk_auto.delta_N, 1)
        self.assertFalse(vtk_auto.enabled)
        self.assertEqual(set(vtk_auto.observables), {"density"})
        self.assertIn(
            f"write to '{key_auto}' every 1 EK steps (disabled)>", repr(vtk_auto))
        key_manual = f"vtk_out/manual_ek_{vtk_suffix}"
        vtk_manual = ek_species.vtk_writers[1]
        self.assertIsInstance(vtk_manual, espressomd.electrokinetics.VTKOutput)
        self.assertEqual(vtk_manual.vtk_uid, key_manual)
        self.assertEqual(vtk_manual.delta_N, 0)
        self.assertEqual(set(vtk_manual.observables), {"density"})
        self.assertIn(f"write to '{key_manual}' on demand>", repr(vtk_manual))
        # check file numbering when resuming VTK write operations
        vtk_root = pathlib.Path("vtk_out") / f"manual_ek_{vtk_suffix}"
        filename = "simulation_step_{}.vtu"
        self.assertTrue((vtk_root / filename.format(0)).exists())
        self.assertFalse((vtk_root / filename.format(1)).exists())
        self.assertFalse((vtk_root / filename.format(2)).exists())
        # check VTK objects are still synchronized with their EK objects
        old_density = ek_species[0, 0, 0].density
        new_density = 1.5 * old_density
        ek_species[0, 0, 0].density = new_density
        vtk_manual.write()
        ek_species[0, 0, 0].density = old_density
        self.assertTrue((vtk_root / filename.format(0)).exists())
        self.assertTrue((vtk_root / filename.format(1)).exists())
        self.assertFalse((vtk_root / filename.format(2)).exists())
        if "espressomd.io.vtk" in sys.modules:
            vtk_reader = espressomd.io.vtk.VTKReader()
            vtk_data = vtk_reader.parse(vtk_root / filename.format(1))
            ek_density = vtk_data["density"]
            self.assertAlmostEqual(
                ek_density[0, 0, 0], new_density, delta=1e-5)

    def test_system_variables(self):
        cell_system_params = system.cell_system.get_state()
        self.assertTrue(cell_system_params['use_verlet_lists'])
        self.assertAlmostEqual(system.cell_system.skin, 0.1, delta=1E-10)
        self.assertAlmostEqual(system.time_step, 0.01, delta=1E-10)
        self.assertAlmostEqual(system.time, 1.5, delta=1E-10)
        self.assertAlmostEqual(system.force_cap, 1e8, delta=1E-10)
        self.assertAlmostEqual(system.min_global_cut, 2.0, delta=1E-10)
        self.assertEqual(system.max_oif_objects, 5)
        np.testing.assert_allclose(np.copy(system.box_l), self.ref_box_l)
        np.testing.assert_array_equal(
            np.copy(system.periodicity), self.ref_periodicity)

    @ut.skipIf('INT.NPT' in modes, 'Lees-Edwards not compatible with NPT')
    def test_lees_edwards(self):
        lebc = system.lees_edwards
        protocol = lebc.protocol
        self.assertEqual(lebc.shear_direction, "x")
        self.assertEqual(lebc.shear_plane_normal, "y")
        self.assertIsInstance(protocol, espressomd.lees_edwards.LinearShear)
        self.assertAlmostEqual(protocol.initial_pos_offset, 0.1, delta=1e-10)
        self.assertAlmostEqual(protocol.time_0, 0.2, delta=1e-10)
        self.assertAlmostEqual(protocol.shear_velocity, 1.2, delta=1e-10)

    def test_particle_properties(self):
        p1, p2, p3, p4, p8 = system.part.by_ids([0, 1, 3, 4, 8])
        np.testing.assert_allclose(np.copy(p1.pos), np.array([1.0, 2.0, 3.0]))
        np.testing.assert_allclose(np.copy(p2.pos), np.array([1.0, 1.0, 2.0]))
        np.testing.assert_allclose(np.copy(p1.f), particle_force0)
        np.testing.assert_allclose(np.copy(p2.f), particle_force1)
        self.assertEqual(p1.type, 0)
        self.assertEqual(p2.type, 0)
        self.assertEqual(p3.type, 1)
        self.assertEqual(p4.type, 1)
        np.testing.assert_allclose(np.copy(p3.v), [0., 0., 0.])
        np.testing.assert_allclose(np.copy(p4.v), [-1., 2., -4.])
        np.testing.assert_allclose(p8.lees_edwards_offset, 0.2)
        np.testing.assert_equal(np.copy(p1.image_box), [0, 0, 0])
        np.testing.assert_equal(np.copy(p8.image_box),
                                np.copy(system.periodicity).astype(int))
        self.assertEqual(p8.lees_edwards_flag, 0)
        if espressomd.has_features('MASS'):
            np.testing.assert_allclose(p3.mass, 1.5)
            np.testing.assert_allclose(p4.mass, 1.)
        if espressomd.has_features('ROTATIONAL_INERTIA'):
            np.testing.assert_allclose(np.copy(p3.rinertia), [2., 3., 4.])
            np.testing.assert_allclose(np.copy(p4.rinertia), [1., 1., 1.])
        if espressomd.has_features('ELECTROSTATICS'):
            np.testing.assert_allclose(p1.q, 1.)
            if espressomd.has_features(['MASS', 'ROTATION']):
                # check Drude particles
                p5 = system.part.by_id(5)
                np.testing.assert_allclose(p2.q, +0.118, atol=1e-3)
                np.testing.assert_allclose(p5.q, -1.118, atol=1e-3)
                np.testing.assert_allclose(p2.mass, 0.4)
                np.testing.assert_allclose(p5.mass, 0.6)
            else:
                np.testing.assert_allclose(p2.q, -1.)
                np.testing.assert_allclose(p2.mass, 1.)
            np.testing.assert_allclose(p3.q, 0.)
            np.testing.assert_allclose(p4.q, 0.)
        if espressomd.has_features('DIPOLES'):
            np.testing.assert_allclose(np.copy(p1.dip), [1.3, 2.1, -6.])
        if espressomd.has_features('ROTATION'):
            np.testing.assert_allclose(np.copy(p3.quat), [1., 2., 3., 4.])
            np.testing.assert_allclose(np.copy(p4.director) * np.sqrt(14.),
                                       [3., 2., 1.])
            np.testing.assert_allclose(np.copy(p3.omega_body), [0., 0., 0.])
            np.testing.assert_allclose(np.copy(p4.omega_body), [0.3, 0.5, 0.7])
            np.testing.assert_equal(np.copy(p3.rotation), [True, False, True])
        if espressomd.has_features('EXTERNAL_FORCES'):
            np.testing.assert_equal(np.copy(p3.fix), [False, True, False])
            np.testing.assert_allclose(np.copy(p3.ext_force), [-0.6, 0.1, 0.2])
        if espressomd.has_features(['EXTERNAL_FORCES', 'ROTATION']):
            np.testing.assert_allclose(np.copy(p3.ext_torque), [0.3, 0.5, 0.7])
        if espressomd.has_features('ROTATIONAL_INERTIA'):
            np.testing.assert_allclose(p3.rinertia, [2., 3., 4.])
        if espressomd.has_features('THERMOSTAT_PER_PARTICLE'):
            gamma = 2.
            if espressomd.has_features('PARTICLE_ANISOTROPY'):
                gamma = np.array([2., 3., 4.])
            np.testing.assert_allclose(p4.gamma, gamma)
            if espressomd.has_features('ROTATION'):
                np.testing.assert_allclose(p3.gamma_rot, 2. * gamma)
        if espressomd.has_features('ENGINE'):
            self.assertEqual(p3.swimming, {"f_swim": 0.03, "mode": "N/A",
                                           "v_swim": 0., "dipole_length": 0.})
        if espressomd.has_features('ENGINE') and has_lb_mode:
            self.assertEqual(p4.swimming, {"v_swim": 0.02, "mode": "puller",
                                           "f_swim": 0., "dipole_length": 1.})
        if espressomd.has_features('LB_ELECTROHYDRODYNAMICS') and has_lb_mode:
            np.testing.assert_allclose(np.copy(p8.mu_E), [-0.1, 0.2, -0.3])
        if espressomd.has_features('VIRTUAL_SITES_RELATIVE'):
            from scipy.spatial.transform import Rotation as R
            q_ind = ([1, 2, 3, 0],)  # convert from scalar-first to scalar-last
            vs_id, vs_dist, vs_quat = p2.vs_relative
            d = p2.pos - p1.pos
            theta = np.arccos(d[2] / np.linalg.norm(d))
            assert abs(theta - 3. * np.pi / 4.) < 1e-8
            q = np.array([0., 0., np.sin(theta / 2.), -np.cos(theta / 2.)])
            r = R.from_quat(p1.quat[q_ind]) * R.from_quat(vs_quat[q_ind])
            self.assertEqual(vs_id, p1.id)
            np.testing.assert_allclose(vs_dist, np.sqrt(2.))
            np.testing.assert_allclose(q[q_ind], r.as_quat(), atol=1e-10)
            np.testing.assert_allclose(np.copy(p2.vs_quat), [1., 0., 0., 0.])

    def test_part_slice(self):
        np.testing.assert_allclose(np.copy(p_slice.id), [4, 1])
        np.testing.assert_allclose(np.copy(p_slice.pos),
                                   np.copy(system.part.by_ids([4, 1]).pos))

    def test_bonded_interactions_serialization(self):
        '''
        Check that particles at the interface between two MPI nodes still
        experience the force from a harmonic bond. The thermostat friction
        is negligible compared to the harmonic bond strength.
        '''
        p3, p4 = system.part.by_ids([3, 4])
        np.testing.assert_allclose(np.copy(p3.pos), system.box_l / 2. - 1.)
        np.testing.assert_allclose(np.copy(p4.pos), system.box_l / 2. + 1.)
        np.testing.assert_allclose(np.copy(p3.f), -np.copy(p4.f), rtol=1e-4)

    @utx.skipIfMissingFeatures('LENNARD_JONES')
    @ut.skipIf('LJ' not in modes, "Skipping test due to missing mode.")
    def test_shape_based_constraints_serialization(self):
        '''
        Check that particles at the interface between two MPI nodes still
        experience the force from a shape-based constraint. The thermostat
        friction is negligible compared to the LJ force.
        '''
        self.assertGreaterEqual(len(system.constraints), 1)
        p3, p4 = system.part.by_ids([3, 4])
        old_force = np.copy(p3.f)
        system.constraints.remove(system.constraints[0])
        old_integrator = system.integrator.integrator
        system.integrator.set_vv()
        system.integrator.run(0, recalc_forces=True)
        system.integrator.integrator = old_integrator
        np.testing.assert_allclose(
            np.copy(p3.f), -np.copy(p4.f), rtol=1e-4)
        self.assertGreater(np.linalg.norm(np.copy(p3.f) - old_force), 1e6)

    @utx.skipIfMissingFeatures(["WALBERLA"])
    @ut.skipIf(not has_lb_mode, "Skipping test due to missing LB feature.")
    @ut.skipIf('THERM.LB' not in modes, 'LB thermostat not in modes')
    def test_thermostat_LB(self):
        thmst = system.thermostat.get_state()[0]
        self.assertEqual(thmst['type'], 'LB')
        # rng_counter_fluid = seed, seed is 0 because kT=0
        self.assertEqual(thmst['rng_counter_fluid'], 0)
        self.assertEqual(thmst['gamma'], 2.0)

    @ut.skipIf('THERM.LANGEVIN' not in modes,
               'Langevin thermostat not in modes')
    def test_thermostat_Langevin(self):
        thmst = system.thermostat.get_state()[0]
        self.assertEqual(thmst['type'], 'LANGEVIN')
        self.assertEqual(thmst['kT'], 1.0)
        self.assertEqual(thmst['seed'], 42)
        self.assertEqual(thmst['counter'], 0)
        self.assertFalse(thmst['act_on_virtual'])
        np.testing.assert_array_equal(thmst['gamma'], 3 * [2.0])
        if espressomd.has_features('ROTATION'):
            np.testing.assert_array_equal(thmst['gamma_rotation'], 3 * [2.0])

    @ut.skipIf('THERM.BD' not in modes,
               'Brownian thermostat not in modes')
    def test_thermostat_Brownian(self):
        thmst = system.thermostat.get_state()[0]
        self.assertEqual(thmst['type'], 'BROWNIAN')
        self.assertEqual(thmst['kT'], 1.0)
        self.assertEqual(thmst['seed'], 42)
        self.assertEqual(thmst['counter'], 0)
        self.assertFalse(thmst['act_on_virtual'])
        np.testing.assert_array_equal(thmst['gamma'], 3 * [2.0])
        if espressomd.has_features('ROTATION'):
            np.testing.assert_array_equal(thmst['gamma_rotation'], 3 * [2.0])

    @utx.skipIfMissingFeatures('DPD')
    @ut.skipIf('THERM.DPD' not in modes, 'DPD thermostat not in modes')
    def test_thermostat_DPD(self):
        thmst = system.thermostat.get_state()[0]
        self.assertEqual(thmst['type'], 'DPD')
        self.assertEqual(thmst['kT'], 1.0)
        self.assertEqual(thmst['seed'], 42)
        self.assertEqual(thmst['counter'], 0)

    @utx.skipIfMissingFeatures('NPT')
    @ut.skipIf('THERM.NPT' not in modes, 'NPT thermostat not in modes')
    def test_thermostat_NPT(self):
        thmst = system.thermostat.get_state()[0]
        self.assertEqual(thmst['type'], 'NPT_ISO')
        self.assertEqual(thmst['seed'], 42)
        self.assertEqual(thmst['counter'], 0)
        self.assertEqual(thmst['gamma0'], 2.0)
        self.assertEqual(thmst['gammav'], 0.1)

    @utx.skipIfMissingFeatures('STOKESIAN_DYNAMICS')
    @ut.skipIf('THERM.SDM' not in modes, 'SDM thermostat not in modes')
    def test_thermostat_SDM(self):
        thmst = system.thermostat.get_state()[0]
        self.assertEqual(thmst['type'], 'SD')
        self.assertEqual(thmst['kT'], 1.0)
        self.assertEqual(thmst['seed'], 42)
        self.assertEqual(thmst['counter'], 0)

    def test_integrator(self):
        params = system.integrator.get_params()
        self.assertAlmostEqual(params['force_cap'], 1e8, delta=1E-10)
        self.assertAlmostEqual(params['time_step'], 0.01, delta=1E-10)
        self.assertAlmostEqual(params['time'], 1.5, delta=1E-10)

    @utx.skipIfMissingFeatures('NPT')
    @ut.skipIf('INT.NPT' not in modes, 'NPT integrator not in modes')
    def test_integrator_NPT(self):
        integ = system.integrator.get_params()
        self.assertIsInstance(
            integ['integrator'], espressomd.integrate.VelocityVerletIsotropicNPT)
        params = integ['integrator'].get_params()
        self.assertEqual(params['ext_pressure'], 2.0)
        self.assertEqual(params['piston'], 0.01)
        self.assertEqual(list(params['direction']), [True, False, False])
        self.assertEqual(params['cubic_box'], False)

    @ut.skipIf('INT.SD' not in modes, 'SD integrator not in modes')
    def test_integrator_SD(self):
        integ = system.integrator.get_params()
        self.assertIsInstance(
            integ['integrator'],
            espressomd.integrate.SteepestDescent)
        params = integ['integrator'].get_params()
        self.assertEqual(params['f_max'], 2.0)
        self.assertEqual(params['gamma'], 0.1)
        self.assertEqual(params['max_displacement'], 0.01)

    @ut.skipIf('INT.NVT' not in modes, 'NVT integrator not in modes')
    def test_integrator_NVT(self):
        integ = system.integrator.get_params()
        self.assertIsInstance(
            integ['integrator'],
            espressomd.integrate.VelocityVerlet)
        params = integ['integrator'].get_params()
        self.assertEqual(params, {})

    @ut.skipIf('INT' in modes, 'VV integrator not the default')
    def test_integrator_VV(self):
        integ = system.integrator.get_params()
        self.assertIsInstance(
            integ['integrator'],
            espressomd.integrate.VelocityVerlet)
        params = integ['integrator'].get_params()
        self.assertEqual(params, {})

    @ut.skipIf('INT.BD' not in modes, 'BD integrator not in modes')
    def test_integrator_BD(self):
        integ = system.integrator.get_params()
        self.assertIsInstance(
            integ['integrator'],
            espressomd.integrate.BrownianDynamics)
        params = integ['integrator'].get_params()
        self.assertEqual(params, {})

    @utx.skipIfMissingFeatures('STOKESIAN_DYNAMICS')
    @ut.skipIf('INT.SDM' not in modes, 'SDM integrator not in modes')
    def test_integrator_SDM(self):
        integ = system.integrator.get_params()
        self.assertIsInstance(
            integ['integrator'],
            espressomd.integrate.StokesianDynamics)
        expected_params = {
            'approximation_method': 'ft', 'radii': {0: 1.5}, 'viscosity': 0.5,
            'lubrication': False, 'pair_mobility': False, 'self_mobility': True}
        params = integ['integrator'].get_params()
        self.assertEqual(params, expected_params)

    @utx.skipIfMissingFeatures('LENNARD_JONES')
    @ut.skipIf('LJ' not in modes, "Skipping test due to missing mode.")
    def test_non_bonded_inter_lj(self):
        self.assertTrue(
            system.non_bonded_inter[0, 0].lennard_jones.call_method("is_registered"))
        params1 = system.non_bonded_inter[0, 0].lennard_jones.get_params()
        params2 = system.non_bonded_inter[3, 0].lennard_jones.get_params()
        reference1 = {'shift': 0.1, 'sigma': 1.3, 'epsilon': 1.2,
                      'cutoff': 2.0, 'offset': 0.0, 'min': 0.0}
        reference2 = {'shift': 0.1, 'sigma': 1.7, 'epsilon': 1.2,
                      'cutoff': 2.0, 'offset': 0.0, 'min': 0.0}
        self.assertEqual(params1, reference1)
        self.assertEqual(params2, reference2)
        self.assertTrue(handle_ia.lennard_jones.call_method("is_registered"))
        self.assertEqual(handle_ia.lennard_jones.get_params(), reference1)

    @utx.skipIfMissingFeatures('DPD')
    def test_non_bonded_inter_dpd(self):
        self.assertEqual(dpd_ia.get_params(), dpd_params)
        self.assertFalse(dpd_ia.call_method("is_registered"))

    def test_bonded_inter(self):
        # check the ObjectHandle was correctly initialized (including MPI)
        bond_ids = system.bonded_inter.call_method('get_bond_ids')
        self.assertEqual(len(bond_ids), len(system.bonded_inter))
        # check bonded interactions
        partcl_1 = system.part.by_id(1)
        reference = {'r_0': 0.0, 'k': 1.0, 'r_cut': 0.0}
        self.assertEqual(partcl_1.bonds[0][0].params, reference)
        self.assertEqual(system.bonded_inter[0].params, reference)
        # all thermalized bonds should be identical
        reference = {**therm_params, 'seed': 3}
        self.assertEqual(partcl_1.bonds[1][0].params, reference)
        self.assertEqual(system.bonded_inter[1].params, reference)
        self.assertEqual(therm_bond2.params, reference)
        # immersed boundary bonds
        self.assertEqual(
            ibm_volcons_bond.params, {'softID': 15, 'kappaV': 0.01})
        self.assertEqual(
            {**ibm_tribend_bond.params, **{'theta0': 0.}},
            {'kb': 2., 'theta0': 0., 'refShape': 'Initial'})
        self.assertEqual(
            ibm_triel_bond.params,
            {'k1': 1.1, 'k2': 1.2, 'maxDist': 1.6, 'elasticLaw': 'NeoHookean'})
        # check new bonds can be added
        if not has_lb_mode:
            new_bond = espressomd.interactions.HarmonicBond(r_0=0.2, k=1.)
            system.bonded_inter.add(new_bond)
            bond_ids = system.bonded_inter.call_method('get_bond_ids')
            self.assertEqual(len(bond_ids), len(system.bonded_inter))

    def test_bond_breakage_specs(self):
        # check the ObjectHandle was correctly initialized (including MPI)
        spec_ids = list(system.bond_breakage.keys())
        self.assertEqual(len(spec_ids), 1)
        cpt_spec = system.bond_breakage[spec_ids[0]]
        self.assertAlmostEqual(
            break_spec.breakage_length,
            cpt_spec.breakage_length,
            delta=1e-10)
        self.assertEqual(break_spec.action_type, cpt_spec.action_type)

    @utx.skipIfMissingFeatures(['ELECTROSTATICS', 'MASS', 'ROTATION'])
    def test_drude_helpers(self):
        drude_type = 10
        core_type = 0
        self.assertIn(drude_type, dh.drude_dict)
        self.assertEqual(dh.drude_dict[drude_type]['alpha'], 1.)
        self.assertEqual(dh.drude_dict[drude_type]['thole_damping'], 2.)
        self.assertEqual(dh.drude_dict[drude_type]['core_type'], core_type)
        self.assertIn(core_type, dh.drude_dict)
        self.assertEqual(dh.drude_dict[core_type]['alpha'], 1.)
        self.assertEqual(dh.drude_dict[core_type]['thole_damping'], 2.)
        self.assertEqual(dh.drude_dict[core_type]['drude_type'], drude_type)
        self.assertEqual(len(dh.drude_dict), 2)
        self.assertEqual(dh.core_type_list, [core_type])
        self.assertEqual(dh.drude_type_list, [drude_type])
        self.assertEqual(dh.core_id_from_drude_id, {5: 1})
        self.assertEqual(dh.drude_id_list, [5])

    @utx.skipIfMissingFeatures(['VIRTUAL_SITES', 'VIRTUAL_SITES_RELATIVE'])
    def test_virtual_sites(self):
        self.assertIsInstance(
            system.virtual_sites,
            espressomd.virtual_sites.VirtualSitesRelative)
        self.assertTrue(system.virtual_sites.have_quaternion)
        self.assertTrue(system.virtual_sites.override_cutoff_check)
        p_real = system.part.by_id(0)
        p_virt = system.part.by_id(1)
        self.assertTrue(p_virt.virtual)
        self.assertFalse(p_real.virtual)
        self.assertEqual(p_real.vs_relative[0], -1)
        self.assertEqual(p_virt.vs_relative[0], p_real.id)
        self.assertEqual(p_real.vs_relative[1], 0.)
        self.assertEqual(p_virt.vs_relative[1], np.sqrt(2.))
        np.testing.assert_allclose(
            p_real.vs_relative[2], [1., 0., 0., 0.], atol=1e-10)

    def test_mean_variance_calculator(self):
        np.testing.assert_array_equal(
            acc_mean_variance.mean(),
            np.array([[1.0, 1.5, 2.0], [1.0, 1.0, 2.0]]))
        np.testing.assert_array_equal(
            acc_mean_variance.variance(),
            np.array([[0., 0.5, 2.], [0., 0., 0.]]))
        np.testing.assert_array_equal(
            system.auto_update_accumulators[0].variance(),
            np.array([[0., 0.5, 2.], [0., 0., 0.]]))

    def test_time_series(self):
        expected = [[[1, 1, 1], [1, 1, 2]], [[1, 2, 3], [1, 1, 2]]]
        np.testing.assert_array_equal(acc_time_series.time_series(), expected)
        np.testing.assert_array_equal(
            system.auto_update_accumulators[1].time_series(),
            expected)

    def test_correlator(self):
        expected = np.zeros((36, 2, 3))
        expected[0:2] = [[[1, 2.5, 5], [1, 1, 4]], [[1, 2, 3], [1, 1, 4]]]
        np.testing.assert_array_equal(acc_correlator.result(), expected)
        np.testing.assert_array_equal(
            system.auto_update_accumulators[2].result(),
            expected)

    @utx.skipIfMissingFeatures('H5MD')
    @utx.skipIfMissingModules("h5py")
    def test_h5md(self):
        # check attributes
        file_path = self.path_cpt_root / "test.h5"
        script_path = pathlib.Path(
            __file__).resolve().parent / "save_checkpoint.py"
        self.assertEqual(h5.fields, ['all'])
        self.assertEqual(h5.script_path, str(script_path))
        self.assertEqual(h5.file_path, str(file_path))

        # write new frame
        h5.write()
        h5.flush()
        h5.close()

        with h5py.File(h5.file_path, 'r') as cur:
            # compare frame #0 against frame #1
            def predicate(cur, key):
                np.testing.assert_allclose(cur[key][0], cur[key][1],
                                           err_msg=f"mismatch for '{key}'")
            for key in ('id', 'position', 'image', 'velocity',
                        'species', 'mass', 'charge', 'force'):
                predicate(cur, f'particles/atoms/{key}/value')
            for key in ('offset', 'direction', 'normal'):
                predicate(cur, f'particles/atoms/lees_edwards/{key}/value')
            predicate(cur, 'particles/atoms/box/edges/value')
            predicate(cur, 'connectivity/atoms/value')

            # check stored physical units
            def predicate(key, attribute):
                self.assertEqual(cur[key].attrs['unit'],
                                 getattr(h5_units, attribute).encode('utf-8'))
            predicate('particles/atoms/id/time', 'time')
            predicate('particles/atoms/lees_edwards/offset/value', 'length')
            predicate('particles/atoms/box/edges/value', 'length')
            predicate('particles/atoms/position/value', 'length')
            predicate('particles/atoms/velocity/value', 'velocity')
            predicate('particles/atoms/force/value', 'force')
            predicate('particles/atoms/charge/value', 'charge')
            predicate('particles/atoms/mass/value', 'mass')

    @utx.skipIfMissingFeatures('DP3M')
    @ut.skipIf('DP3M.CPU' not in modes,
               "Skipping test due to missing combination.")
    def test_dp3m(self):
        actor = self.get_active_actor_of_type(
            espressomd.magnetostatics.DipolarP3M)
        state = actor.get_params()
        reference = {'prefactor': 1.0, 'accuracy': 0.01, 'mesh': 3 * [8],
                     'cao': 1, 'alpha': 12.0, 'r_cut': 2.4, 'tune': False,
                     'mesh_off': [0.5, 0.5, 0.5], 'epsilon': 2.0, 'timings': 15}
        for key in reference:
            self.assertIn(key, state)
            np.testing.assert_almost_equal(state[key], reference[key],
                                           err_msg=f'for parameter {key}')

    @utx.skipIfMissingFeatures('P3M')
    @ut.skipIf(not has_p3m_mode, "Skipping test due to missing combination.")
    def test_p3m(self):
        actor = self.get_active_actor_of_type(
            espressomd.electrostatics._P3MBase)
        state = actor.get_params()
        reference = {'prefactor': 1.0, 'accuracy': 0.1, 'mesh': 3 * [10],
                     'cao': 1, 'alpha': 1.0, 'r_cut': 1.0, 'tune': False,
                     'timings': 15, 'check_neutrality': True,
                     'check_complex_residuals': False,
                     'charge_neutrality_tolerance': 1e-12}
        for key in reference:
            self.assertIn(key, state)
            np.testing.assert_almost_equal(state[key], reference[key],
                                           err_msg=f'for parameter {key}')

    @utx.skipIfMissingFeatures('P3M')
    @ut.skipIf('ELC' not in modes, "Skipping test due to missing combination.")
    def test_elc(self):
        actor = self.get_active_actor_of_type(espressomd.electrostatics.ELC)
        elc_state = actor.get_params()
        p3m_state = elc_state['actor'].get_params()
        p3m_reference = {'prefactor': 1.0, 'accuracy': 0.1, 'mesh': 3 * [10],
                         'cao': 1, 'alpha': 1.0, 'r_cut': 1.0, 'tune': False,
                         'timings': 15, 'check_neutrality': True,
                         'check_complex_residuals': False,
                         'charge_neutrality_tolerance': 7e-12}
        elc_reference = {'gap_size': 6.0, 'maxPWerror': 0.1,
                         'delta_mid_top': 0.9, 'delta_mid_bot': 0.1,
                         'check_neutrality': True,
                         'charge_neutrality_tolerance': 5e-12}
        for key in elc_reference:
            self.assertIn(key, elc_state)
            np.testing.assert_almost_equal(elc_state[key], elc_reference[key],
                                           err_msg=f'for parameter {key}')
        for key in p3m_reference:
            self.assertIn(key, p3m_state)
            np.testing.assert_almost_equal(p3m_state[key], p3m_reference[key],
                                           err_msg=f'for parameter {key}')

    @utx.skipIfMissingFeatures(["SCAFACOS"])
    @utx.skipIfMissingScafacosMethod("p3m")
    @ut.skipIf('SCAFACOS' not in modes, "Missing combination.")
    def test_scafacos_coulomb(self):
        actor = self.get_active_actor_of_type(
            espressomd.electrostatics.Scafacos)
        state = actor.get_params()
        reference = {'prefactor': 0.5, 'method_name': 'p3m',
                     'method_params': {
                         'p3m_cao': 7,
                         'p3m_r_cut': 1.0,
                         'p3m_grid': 64,
                         'p3m_alpha': 2.084652}}
        for key in reference:
            self.assertEqual(state[key], reference[key], msg=f'for {key}')

    @utx.skipIfMissingFeatures(["SCAFACOS_DIPOLES"])
    @utx.skipIfMissingScafacosMethod("p2nfft")
    @ut.skipIf('SCAFACOS' not in modes, "Missing combination.")
    def test_scafacos_dipoles(self):
        actor = self.get_active_actor_of_type(
            espressomd.magnetostatics.Scafacos)
        state = actor.get_params()
        reference = {'prefactor': 1.2, 'method_name': 'p2nfft',
                     'method_params': {
                         "p2nfft_verbose_tuning": 0,
                         "pnfft_N": [32, 32, 32],
                         "pnfft_n": [32, 32, 32],
                         "pnfft_window_name": "bspline",
                         "pnfft_m": 4,
                         "p2nfft_ignore_tolerance": 1,
                         "pnfft_diff_ik": 0,
                         "p2nfft_r_cut": 11,
                         "p2nfft_alpha": 0.37}}
        for key in reference:
            self.assertIn(key, state)
            self.assertEqual(state[key], reference[key], msg=f'for {key}')

    def test_comfixed(self):
        self.assertEqual(list(system.comfixed.types), [0, 2])

    @utx.skipIfMissingFeatures('COLLISION_DETECTION')
    def test_collision_detection(self):
        coldet = system.collision_detection
        self.assertEqual(coldet.mode, "bind_centers")
        self.assertAlmostEqual(coldet.distance, 0.11, delta=1E-9)
        self.assertEqual(coldet.bond_centers, system.bonded_inter[0])

    @utx.skipIfMissingFeatures('EXCLUSIONS')
    def test_exclusions(self):
        self.assertEqual(list(system.part.by_id(0).exclusions), [2])
        self.assertEqual(list(system.part.by_id(1).exclusions), [2])
        self.assertEqual(list(system.part.by_id(2).exclusions), [0, 1])

    def test_constraints(self):
        n_contraints = 8
        if espressomd.has_features("ELECTROSTATICS"):
            n_contraints += 1
        self.assertEqual(len(system.constraints), n_contraints)

        c = system.constraints
        ref_shape = self.ref_box_l.astype(int) + 2

        self.assertIsInstance(c[0].shape, espressomd.shapes.Sphere)
        self.assertAlmostEqual(c[0].shape.radius, 0.1, delta=1E-10)
        self.assertEqual(c[0].particle_type, 7)

        self.assertIsInstance(c[1].shape, espressomd.shapes.Wall)
        np.testing.assert_allclose(np.copy(c[1].shape.normal),
                                   [1. / np.sqrt(3)] * 3)

        self.assertIsInstance(c[2], espressomd.constraints.Gravity)
        np.testing.assert_allclose(np.copy(c[2].g), [1., 2., 3.])

        self.assertIsInstance(
            c[3], espressomd.constraints.HomogeneousMagneticField)
        np.testing.assert_allclose(np.copy(c[3].H), [1., 2., 3.])

        self.assertIsInstance(
            c[4], espressomd.constraints.HomogeneousFlowField)
        np.testing.assert_allclose(np.copy(c[4].u), [1., 2., 3.])
        self.assertAlmostEqual(c[4].gamma, 2.3, delta=1E-10)

        self.assertIsInstance(c[5], espressomd.constraints.PotentialField)
        self.assertEqual(c[5].field.shape, tuple(list(ref_shape) + [1]))
        self.assertAlmostEqual(c[5].default_scale, 1.6, delta=1E-10)
        self.assertAlmostEqual(c[5].particle_scales[5], 6.0, delta=1E-10)
        np.testing.assert_allclose(np.copy(c[5].origin), [-0.5, -0.5, -0.5])
        np.testing.assert_allclose(np.copy(c[5].grid_spacing), np.ones(3))
        ref_pot = espressomd.constraints.PotentialField(
            field=pot_field_data, grid_spacing=np.ones(3), default_scale=1.6)
        np.testing.assert_allclose(np.copy(c[5].field), np.copy(ref_pot.field),
                                   atol=1e-10)

        self.assertIsInstance(c[6], espressomd.constraints.ForceField)
        self.assertEqual(c[6].field.shape, tuple(list(ref_shape) + [3]))
        self.assertAlmostEqual(c[6].default_scale, 1.4, delta=1E-10)
        np.testing.assert_allclose(np.copy(c[6].origin), [-0.5, -0.5, -0.5])
        np.testing.assert_allclose(np.copy(c[6].grid_spacing), np.ones(3))
        ref_vec = espressomd.constraints.ForceField(
            field=vec_field_data, grid_spacing=np.ones(3), default_scale=1.4)
        np.testing.assert_allclose(np.copy(c[6].field), np.copy(ref_vec.field),
                                   atol=1e-10)

        union = c[7].shape
        self.assertIsInstance(union, espressomd.shapes.Union)
        self.assertEqual(c[7].particle_type, 2)
        self.assertEqual(len(union), 2)
        wall1, wall2 = union.call_method('get_elements')
        self.assertIsInstance(wall1, espressomd.shapes.Wall)
        self.assertIsInstance(wall2, espressomd.shapes.Wall)
        np.testing.assert_allclose(np.copy(wall1.normal),
                                   [1., 0., 0.], atol=1e-10)
        np.testing.assert_allclose(np.copy(wall2.normal),
                                   [0., 1., 0.], atol=1e-10)
        np.testing.assert_allclose(wall1.dist, 0.5, atol=1e-10)
        np.testing.assert_allclose(wall2.dist, 1.5, atol=1e-10)

        if espressomd.has_features("ELECTROSTATICS"):
            wave = c[n_contraints - 1]
            self.assertIsInstance(
                wave, espressomd.constraints.ElectricPlaneWave)
            np.testing.assert_allclose(np.copy(wave.E0), [1., -2., 3.])
            np.testing.assert_allclose(np.copy(wave.k), [-.1, .2, .3])
            self.assertAlmostEqual(wave.omega, 5., delta=1E-10)
            self.assertAlmostEqual(wave.phi, 1.4, delta=1E-10)

    @utx.skipIfMissingFeatures("WCA")
    @ut.skipIf(has_lb_mode, "LB not supported")
    @ut.skipIf("INT.SDM" in modes, "Stokesian integrator not supported")
    @ut.skipIf("INT.BD" in modes, "Brownian integrator not supported")
    @ut.skipIf("INT.SD" in modes, "Steepest descent not supported")
    def test_union(self):
        # the union shape is an object list, and should be properly
        # deserialized on all MPI ranks
        system.non_bonded_inter[2, 6].wca.set_params(epsilon=1., sigma=1.)
        p1 = system.part.add(pos=[1., 1.6, 0.], type=6)
        p2 = system.part.add(pos=[system.box_l[0] - 1., 1.6, 0.], type=6)
        system.integrator.run(0, recalc_forces=True)
        np.testing.assert_allclose(p1.f, [0., 1e8, 0.], atol=1e-3)
        np.testing.assert_allclose(p2.f, [0., 1e8, 0.], atol=1e-3)
        p1.remove()
        p2.remove()
        system.non_bonded_inter[2, 6].reset()


if __name__ == '__main__':
    config.bind_test_class(CheckpointTest)
    ut.main()
