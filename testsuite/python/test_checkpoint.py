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
    n_nodes = system.cell_system.get_state()["n_nodes"]

    @classmethod
    def setUpClass(cls):
        cls.ref_box_l = np.array([12.0, 14.0, 16.0])
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

    @utx.skipIfMissingFeatures('WALBERLA')
    @ut.skipIf(not has_lb_mode, "Skipping test due to missing LB feature.")
    def test_lb_fluid(self):
        lbf = self.get_active_actor_of_type(espressomd.lb.LBFluidWalberla)
        cpt_mode = 0 if 'LB.ASCII' in modes else 1
        cpt_root = pathlib.Path(self.checkpoint.checkpoint_dir)
        cpt_path = str(cpt_root / "lb") + "{}.cpt"

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
        precision = 9 if "LB.WALBERLA" in modes else 5
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
            'agrid': 1.0,
            'viscosity': 1.3,
            'density': 1.5,
            'tau': 0.01}
        for key in reference:
            self.assertIn(key, state)
            self.assertAlmostEqual(reference[key], state[key], delta=1E-7)
        self.assertTrue(lbf.is_active)
        self.assertFalse(lbf.is_single_precision)

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

    @utx.skipIfMissingFeatures('WALBERLA')
    @ut.skipIf(not has_lb_mode, "Skipping test due to missing LB feature.")
    def test_lb_vtk(self):
        vtk_suffix = config.test_name
        vtk_registry = espressomd.lb._vtk_registry
        key_auto = f'vtk_out/auto_lb_{vtk_suffix}'
        self.assertIn(key_auto, vtk_registry.map)
        obj = vtk_registry.map[key_auto]
        self.assertIsInstance(obj, espressomd.lb.VTKOutput)
        self.assertEqual(obj.vtk_uid, key_auto)
        self.assertEqual(obj.delta_N, 1)
        self.assertFalse(obj.enabled)
        self.assertEqual(set(obj.observables), {'density', 'velocity_vector'})
        self.assertIn(
            f"write to '{key_auto}' every 1 LB steps (disabled)>", repr(obj))
        key_manual = f'vtk_out/manual_lb_{vtk_suffix}'
        self.assertIn(key_manual, vtk_registry.map)
        obj = vtk_registry.map[key_manual]
        self.assertIsInstance(obj, espressomd.lb.VTKOutput)
        self.assertEqual(obj.vtk_uid, key_manual)
        self.assertEqual(obj.delta_N, 0)
        self.assertEqual(set(obj.observables), {'density'})
        self.assertIn(f"write to '{key_manual}' on demand>", repr(obj))
        # check file numbering when resuming VTK write operations
        vtk_root = pathlib.Path("vtk_out") / f"manual_lb_{vtk_suffix}"
        filename = "simulation_step_{}.vtu"
        vtk_manual = espressomd.lb._vtk_registry.map[key_manual]
        self.assertTrue((vtk_root / filename.format(0)).exists())
        self.assertFalse((vtk_root / filename.format(1)).exists())
        self.assertFalse((vtk_root / filename.format(2)).exists())
        vtk_manual.write()
        self.assertTrue((vtk_root / filename.format(0)).exists())
        self.assertTrue((vtk_root / filename.format(1)).exists())
        self.assertFalse((vtk_root / filename.format(2)).exists())

    # TODO walberla
#    @utx.skipIfMissingFeatures('WALBERLA')
#    @ut.skipIf(not has_lb_mode, "Skipping test due to missing LB feature.")
#    def test_ek_vtk(self):
#        vtk_suffix = config.test_name
#        vtk_registry = espressomd.EKSpecies._ek_vtk_registry
#        key_auto = f'vtk_out/auto_ek_{vtk_suffix}'
#        self.assertIn(key_auto, vtk_registry.map)
#        obj = vtk_registry.map[key_auto]
#        self.assertIsInstance(obj, espressomd.lb.VTKOutput)
#        self.assertEqual(obj.vtk_uid, key_auto)
#        self.assertEqual(obj.delta_N, 1)
#        self.assertFalse(obj.enabled)
#        self.assertEqual(set(obj.observables), {'density'})
#        self.assertIn(
#            f"write to '{key_auto}' every 1 LB steps (disabled)>", repr(obj))
#        key_manual = f'vtk_out/manual_ek_{vtk_suffix}'
#        self.assertIn(key_manual, vtk_registry.map)
#        obj = vtk_registry.map[key_manual]
#        self.assertIsInstance(obj, espressomd.lb.VTKOutput)
#        self.assertEqual(obj.vtk_uid, key_manual)
#        self.assertEqual(obj.delta_N, 0)
#        self.assertEqual(set(obj.observables), {'density'})
#        self.assertIn(f"write to '{key_manual}' on demand>", repr(obj))
#        # check file numbering when resuming VTK write operations
#        vtk_root = pathlib.Path("vtk_out") / f"manual_ek_{vtk_suffix}"
#        filename = "simulation_step_{}.vtu"
#        vtk_manual = espressomd.lb._vtk_registry.map[key_manual]
#        self.assertTrue((vtk_root / filename.format(0)).exists())
#        self.assertFalse((vtk_root / filename.format(1)).exists())
#        self.assertFalse((vtk_root / filename.format(2)).exists())
#        vtk_manual.write()
#        self.assertTrue((vtk_root / filename.format(0)).exists())
#        self.assertTrue((vtk_root / filename.format(1)).exists())
#        self.assertFalse((vtk_root / filename.format(2)).exists())

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

    def test_part(self):
        p1, p2 = system.part.by_ids([0, 1])
        np.testing.assert_allclose(np.copy(p1.pos), np.array([1.0, 2.0, 3.0]))
        np.testing.assert_allclose(np.copy(p2.pos), np.array([1.0, 1.0, 2.0]))
        np.testing.assert_allclose(np.copy(p1.f), particle_force0)
        np.testing.assert_allclose(np.copy(p2.f), particle_force1)

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
        system.integrator.run(0, recalc_forces=True)
        np.testing.assert_allclose(
            np.copy(p3.f), -np.copy(p4.f), rtol=1e-4)
        self.assertGreater(np.linalg.norm(np.copy(p3.f) - old_force), 1e6)

    @utx.skipIfMissingFeatures('WALBERLA')
    @ut.skipIf(not has_lb_mode, "Skipping test due to missing LB feature.")
    @ut.skipIf('THERM.LB' not in modes, 'LB thermostat not in modes')
    def test_thermostat_LB(self):
        thmst = system.thermostat.get_state()[0]
        # TODO WALBERLA
#        if 'LB.GPU' in modes and not espressomd.gpu_available():
#            self.assertEqual(thmst['type'], 'OFF')
#        else:
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
        if 'THERM.LB' not in modes:
            # all thermalized bonds should be identical
            reference = {**thermalized_bond_params, 'seed': 3}
            self.assertEqual(partcl_1.bonds[1][0].params, reference)
            self.assertEqual(system.bonded_inter[1].params, reference)
            self.assertEqual(thermalized_bond2.params, reference)
        # immersed boundary bonds
        self.assertEqual(
            ibm_volcons_bond.params, {'softID': 15, 'kappaV': 0.01})
        if 'DP3M.CPU' not in modes:
            self.assertEqual(
                ibm_tribend_bond.params,
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

    @ut.skipIf('THERM.LB' in modes, 'LB thermostat in modes')
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
        self.assertTrue(system.part.by_id(1).virtual)
        self.assertIsInstance(
            system.virtual_sites,
            espressomd.virtual_sites.VirtualSitesRelative)

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
            # note: cannot compare forces since they are NaN in frame #1
            for key in ('position', 'image', 'velocity',
                        'id', 'species', 'mass', 'charge'):
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
        n_contraints = 7
        if self.n_nodes == 1:
            n_contraints += 1
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

        if self.n_nodes == 1:
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


if __name__ == '__main__':
    config.bind_test_class(CheckpointTest)
    ut.main()
