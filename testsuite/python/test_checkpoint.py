# Copyright (C) 2010-2019 The ESPResSo project
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

# pylint: disable=undefined-variable
import unittest as ut
import unittest_decorators as utx
import numpy as np

import espressomd
import espressomd.checkpointing
import espressomd.electrostatics
import espressomd.magnetostatics
import espressomd.scafacos
import espressomd.virtual_sites
import espressomd.integrate
import espressomd.shapes
import espressomd.constraints

modes = {x for mode in set("@TEST_COMBINATION@".upper().split('-'))
         for x in [mode, mode.split('.')[0]]}

LB = ('LB.CPU' in modes or 'LB.GPU' in modes and espressomd.gpu_available())


class CheckpointTest(ut.TestCase):

    checkpoint = espressomd.checkpointing.Checkpoint(
        checkpoint_id="mycheckpoint_@TEST_COMBINATION@_@TEST_BINARY@".replace(
            '.', '__'),
        checkpoint_path="@CMAKE_CURRENT_BINARY_DIR@")
    checkpoint.load(0)
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

    @ut.skipIf(not LB, "Skipping test due to missing mode.")
    def test_lb_fluid(self):
        '''
        Check serialization of the LB fluid. The checkpoint file only stores
        population information, therefore calling ``lbf.load_checkpoint()``
        erases all LBBoundaries information but doesn't remove the objects
        contained in ``system.lbboundaries`. This test method is named such
        that it is executed after ``self.test_lb_boundaries()``.
        '''
        lbf = self.get_active_actor_of_type(
            espressomd.lb.HydrodynamicInteraction)
        cpt_mode = int("@TEST_BINARY@")
        cpt_path = self.checkpoint.checkpoint_dir + "/lb{}.cpt"

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
        with self.assertRaisesRegex(RuntimeError, 'could not open file'):
            lbf.load_checkpoint(cpt_path.format("-unknown"), cpt_mode)

        # load the valid LB checkpoint file
        lbf.load_checkpoint(cpt_path.format(""), cpt_mode)
        precision = 9 if "LB.CPU" in modes else 5
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
        state = lbf.get_params()
        reference = {'agrid': 0.5, 'visc': 1.3, 'dens': 1.5, 'tau': 0.01,
                     'gamma_odd': 0.2, 'gamma_even': 0.3}
        for key in reference:
            self.assertIn(key, state)
            self.assertAlmostEqual(reference[key], state[key], delta=1E-7)

    def test_system_variables(self):
        cell_system_params = system.cell_system.get_state()
        self.assertTrue(cell_system_params['use_verlet_list'])
        self.assertAlmostEqual(system.cell_system.skin, 0.1, delta=1E-10)
        self.assertAlmostEqual(system.time_step, 0.01, delta=1E-10)
        self.assertAlmostEqual(system.time, 1.5, delta=1E-10)
        self.assertAlmostEqual(system.force_cap, 1e8, delta=1E-10)
        self.assertAlmostEqual(system.min_global_cut, 2.0, delta=1E-10)
        self.assertEqual(system.max_oif_objects, 5)
        np.testing.assert_allclose(np.copy(system.box_l), self.ref_box_l)
        np.testing.assert_array_equal(
            np.copy(system.periodicity), self.ref_periodicity)

    def test_part(self):
        p1, p2 = system.part[0:2]
        np.testing.assert_allclose(np.copy(p1.pos), np.array([1.0, 2.0, 3.0]))
        np.testing.assert_allclose(np.copy(p2.pos), np.array([1.0, 1.0, 2.0]))
        np.testing.assert_allclose(np.copy(p1.f), particle_force0)
        np.testing.assert_allclose(np.copy(p2.f), particle_force1)

    def test_bonded_interactions_serialization(self):
        '''
        Check that particles at the interface between two MPI nodes still
        experience the force from a harmonic bond. The thermostat friction
        is negligible compared to the harmonic bond strength.
        '''
        p3, p4 = system.part[3:5]
        np.testing.assert_allclose(np.copy(p3.pos), system.box_l / 2. - 1.)
        np.testing.assert_allclose(np.copy(p4.pos), system.box_l / 2. + 1.)
        np.testing.assert_allclose(np.copy(p3.f), -np.copy(p4.f), rtol=1e-4)

    @utx.skipIfMissingFeatures('LENNARD_JONES')
    @ut.skipIf('LJ' not in modes, "Skipping test due to missing mode.")
    @ut.skipIf(n_nodes > 1, "only runs for 1 MPI rank")
    def test_shape_based_constraints_serialization(self):
        '''
        Check that particles at the interface between two MPI nodes still
        experience the force from a shape-based constraint. The thermostat
        friction is negligible compared to the LJ force.
        '''
        p3, p4 = system.part[3:5]
        old_force = np.copy(p3.f)
        system.constraints.remove(system.constraints[0])
        system.integrator.run(0, recalc_forces=True)
        np.testing.assert_allclose(
            np.copy(p3.f), -np.copy(p4.f), rtol=1e-4)
        self.assertGreater(np.linalg.norm(np.copy(p3.f) - old_force), 1e6)

    @ut.skipIf('THERM.LB' not in modes, 'LB thermostat not in modes')
    def test_thermostat_LB(self):
        thmst = system.thermostat.get_state()[0]
        if 'LB.GPU' in modes and not espressomd.gpu_available():
            self.assertEqual(thmst['type'], 'OFF')
        else:
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
        params = system.integrator.get_state()
        self.assertAlmostEqual(params['force_cap'], 1e8, delta=1E-10)
        self.assertAlmostEqual(params['time_step'], 0.01, delta=1E-10)
        self.assertAlmostEqual(params['time'], 1.5, delta=1E-10)

    @utx.skipIfMissingFeatures('NPT')
    @ut.skipIf('INT.NPT' not in modes, 'NPT integrator not in modes')
    def test_integrator_NPT(self):
        integ = system.integrator.get_state()
        self.assertIsInstance(
            integ['integrator'], espressomd.integrate.VelocityVerletIsotropicNPT)
        params = integ['integrator'].get_params()
        self.assertEqual(params['ext_pressure'], 2.0)
        self.assertEqual(params['piston'], 0.01)
        self.assertEqual(params['direction'], [1, 0, 0])
        self.assertEqual(params['cubic_box'], False)

    @ut.skipIf('INT.SD' not in modes, 'SD integrator not in modes')
    def test_integrator_SD(self):
        integ = system.integrator.get_state()
        self.assertIsInstance(
            integ['integrator'],
            espressomd.integrate.SteepestDescent)
        params = integ['integrator'].get_params()
        self.assertEqual(params['f_max'], 2.0)
        self.assertEqual(params['gamma'], 0.1)
        self.assertEqual(params['max_displacement'], 0.01)

    @ut.skipIf('INT.NVT' not in modes, 'NVT integrator not in modes')
    def test_integrator_NVT(self):
        integ = system.integrator.get_state()
        self.assertIsInstance(
            integ['integrator'],
            espressomd.integrate.VelocityVerlet)
        params = integ['integrator'].get_params()
        self.assertEqual(params, {})

    @ut.skipIf('INT' in modes, 'VV integrator not the default')
    def test_integrator_VV(self):
        integ = system.integrator.get_state()
        self.assertIsInstance(
            integ['integrator'],
            espressomd.integrate.VelocityVerlet)
        params = integ['integrator'].get_params()
        self.assertEqual(params, {})

    @ut.skipIf('INT.BD' not in modes, 'BD integrator not in modes')
    def test_integrator_BD(self):
        integ = system.integrator.get_state()
        self.assertIsInstance(
            integ['integrator'],
            espressomd.integrate.BrownianDynamics)
        params = integ['integrator'].get_params()
        self.assertEqual(params, {})

    @utx.skipIfMissingFeatures('STOKESIAN_DYNAMICS')
    @ut.skipIf('INT.SDM' not in modes, 'SDM integrator not in modes')
    def test_integrator_SDM(self):
        integ = system.integrator.get_state()
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
    def test_non_bonded_inter(self):
        params1 = system.non_bonded_inter[0, 0].lennard_jones.get_params()
        params2 = system.non_bonded_inter[3, 0].lennard_jones.get_params()
        reference1 = {'shift': 0.1, 'sigma': 1.3, 'epsilon': 1.2,
                      'cutoff': 2.0, 'offset': 0.0, 'min': 0.0}
        reference2 = {'shift': 0.1, 'sigma': 1.7, 'epsilon': 1.2,
                      'cutoff': 2.0, 'offset': 0.0, 'min': 0.0}
        self.assertEqual(params1, reference1)
        self.assertEqual(params2, reference2)

    def test_bonded_inter(self):
        state = system.part[1].bonds[0][0].params
        reference = {'r_0': 0.0, 'k': 1.0, 'r_cut': 0.0}
        self.assertEqual(state, reference)
        state = system.part[1].bonds[0][0].params
        self.assertEqual(state, reference)
        if 'THERM.LB' not in modes:
            state = system.part[1].bonds[1][0].params
            reference = {'temp_com': 0., 'gamma_com': 0., 'temp_distance': 0.2,
                         'gamma_distance': 0.5, 'r_cut': 2.0, 'seed': 51}
            self.assertEqual(state, reference)
            state = system.part[1].bonds[1][0].params
            self.assertEqual(state, reference)
        # immersed boundary bonds
        self.assertEqual(
            ibm_volcons_bond.params, {'softID': 15, 'kappaV': 0.01})
        if 'DP3M.CPU' not in modes:
            self.assertEqual(ibm_tribend_bond.params, {'kb': 2., 'theta0': 0.})
        self.assertEqual(
            ibm_triel_bond.params,
            {'k1': 1.1, 'k2': 1.2, 'maxDist': 1.6, 'elasticLaw': 'NeoHookean'})

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
        self.assertTrue(system.part[1].virtual)
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
        if self.n_nodes == 1:
            np.testing.assert_array_equal(
                system.auto_update_accumulators[0].variance(),
                np.array([[0., 0.5, 2.], [0., 0., 0.]]))

    def test_time_series(self):
        expected = [[[1, 1, 1], [1, 1, 2]], [[1, 2, 3], [1, 1, 2]]]
        np.testing.assert_array_equal(acc_time_series.time_series(), expected)
        if self.n_nodes == 1:
            np.testing.assert_array_equal(
                system.auto_update_accumulators[1].time_series(),
                expected)

    def test_correlator(self):
        expected = np.zeros((36, 2, 3))
        expected[0:2] = [[[1, 2.5, 5], [1, 1, 4]], [[1, 2, 3], [1, 1, 4]]]
        np.testing.assert_array_equal(acc_correlator.result(), expected)
        if self.n_nodes == 1:
            np.testing.assert_array_equal(
                system.auto_update_accumulators[2].result(),
                expected)

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
    @ut.skipIf('P3M.CPU' not in modes,
               "Skipping test due to missing combination.")
    def test_p3m(self):
        actor = self.get_active_actor_of_type(espressomd.electrostatics.P3M)
        state = actor.get_params()
        reference = {'prefactor': 1.0, 'accuracy': 0.1, 'mesh': 3 * [10],
                     'cao': 1, 'alpha': 1.0, 'r_cut': 1.0, 'tune': False,
                     'timings': 15}
        for key in reference:
            self.assertIn(key, state)
            np.testing.assert_almost_equal(state[key], reference[key],
                                           err_msg=f'for parameter {key}')

    @utx.skipIfMissingFeatures('P3M')
    @ut.skipIf('P3M.ELC' not in modes,
               "Skipping test due to missing combination.")
    def test_elc(self):
        actor = self.get_active_actor_of_type(espressomd.electrostatics.ELC)
        elc_state = actor.get_params()
        p3m_state = elc_state['p3m_actor'].get_params()
        p3m_reference = {'prefactor': 1.0, 'accuracy': 0.1, 'mesh': 3 * [10],
                         'cao': 1, 'alpha': 1.0, 'r_cut': 1.0, 'tune': False,
                         'timings': 15}
        elc_reference = {'gap_size': 6.0, 'maxPWerror': 0.1,
                         'delta_mid_top': 0.9, 'delta_mid_bot': 0.1}
        for key in elc_reference:
            self.assertIn(key, elc_state)
            np.testing.assert_almost_equal(elc_state[key], elc_reference[key],
                                           err_msg=f'for parameter {key}')
        for key in p3m_reference:
            self.assertIn(key, p3m_state)
            np.testing.assert_almost_equal(p3m_state[key], p3m_reference[key],
                                           err_msg=f'for parameter {key}')

    @ut.skipIf(not espressomd.has_features('SCAFACOS') or
               'SCAFACOS' not in modes or
               'p3m' not in espressomd.scafacos.available_methods(),
               "Skipping test due to missing combination or p3m method.")
    def test_scafacos(self):
        actor = self.get_active_actor_of_type(
            espressomd.electrostatics.Scafacos)
        state = actor.get_params()
        reference = {'prefactor': 0.5, 'method_name': 'p3m',
                     'method_params': {
                         'p3m_cao': '7',
                         'p3m_r_cut': '1.0',
                         'p3m_grid': '64',
                         'p3m_alpha': '2.084652'}}
        for key in reference:
            self.assertEqual(state[key], reference[key], msg=f'for {key}')

    @ut.skipIf(not espressomd.has_features('SCAFACOS_DIPOLES') or
               'SCAFACOS' not in modes or
               'p2nfft' not in espressomd.scafacos.available_methods(),
               "Skipping test due to missing combination or p2nfft method.")
    def test_scafacos_dipoles(self):
        actor = self.get_active_actor_of_type(
            espressomd.magnetostatics.Scafacos)
        state = actor.get_params()
        reference = {'prefactor': 1.2, 'method_name': 'p2nfft',
                     'method_params': {
                         "p2nfft_verbose_tuning": "0",
                         "pnfft_N": "32,32,32",
                         "pnfft_n": "32,32,32",
                         "pnfft_window_name": "bspline",
                         "pnfft_m": "4",
                         "p2nfft_ignore_tolerance": "1",
                         "pnfft_diff_ik": "0",
                         "p2nfft_r_cut": "11",
                         "p2nfft_alpha": "0.37"}}
        for key in reference:
            self.assertIn(key, state)
            self.assertEqual(state[key], reference[key], msg=f'for {key}')

    @utx.skipIfMissingFeatures('COLLISION_DETECTION')
    def test_collision_detection(self):
        coldet = system.collision_detection
        self.assertEqual(coldet.mode, "bind_centers")
        self.assertAlmostEqual(coldet.distance, 0.11, delta=1E-9)
        self.assertEqual(coldet.bond_centers, system.bonded_inter[0])

    @utx.skipIfMissingFeatures('EXCLUSIONS')
    def test_exclusions(self):
        self.assertEqual(list(system.part[0].exclusions), [2])
        self.assertEqual(list(system.part[1].exclusions), [2])
        self.assertEqual(list(system.part[2].exclusions), [0, 1])

    @ut.skipIf(not LB or not (espressomd.has_features("LB_BOUNDARIES")
                              or espressomd.has_features("LB_BOUNDARIES_GPU")), "Missing features")
    @ut.skipIf(n_nodes > 1, "only runs for 1 MPI rank")
    def test_lb_boundaries(self):
        # check boundary objects
        self.assertEqual(len(system.lbboundaries), 2)
        np.testing.assert_allclose(
            np.copy(system.lbboundaries[0].velocity), [1e-4, 1e-4, 0])
        np.testing.assert_allclose(
            np.copy(system.lbboundaries[1].velocity), [0, 0, 0])
        self.assertIsInstance(
            system.lbboundaries[0].shape, espressomd.shapes.Wall)
        self.assertIsInstance(
            system.lbboundaries[1].shape, espressomd.shapes.Wall)
        # check boundary flag
        lbf = self.get_active_actor_of_type(
            espressomd.lb.HydrodynamicInteraction)
        np.testing.assert_equal(np.copy(lbf[0, :, :].boundary.astype(int)), 1)
        np.testing.assert_equal(np.copy(lbf[-1, :, :].boundary.astype(int)), 2)
        np.testing.assert_equal(
            np.copy(lbf[1:-1, :, :].boundary.astype(int)), 0)
        # remove boundaries
        system.lbboundaries.clear()
        self.assertEqual(len(system.lbboundaries), 0)
        np.testing.assert_equal(np.copy(lbf[:, :, :].boundary.astype(int)), 0)

    @ut.skipIf(n_nodes > 1, "only runs for 1 MPI rank")
    def test_constraints(self):
        self.assertEqual(len(system.constraints),
                         8 - int(not espressomd.has_features("ELECTROSTATICS")))
        c = system.constraints
        ref_shape = self.ref_box_l.astype(int) + 2

        self.assertIsInstance(c[0].shape, espressomd.shapes.Sphere)
        self.assertAlmostEqual(c[0].shape.radius, 0.1, delta=1E-10)
        self.assertEqual(c[0].particle_type, 17)

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

        if espressomd.has_features("ELECTROSTATICS"):
            self.assertIsInstance(
                c[7], espressomd.constraints.ElectricPlaneWave)
            np.testing.assert_allclose(np.copy(c[7].E0), [1., -2., 3.])
            np.testing.assert_allclose(np.copy(c[7].k), [-.1, .2, .3])
            self.assertAlmostEqual(c[7].omega, 5., delta=1E-10)
            self.assertAlmostEqual(c[7].phi, 1.4, delta=1E-10)


if __name__ == '__main__':
    ut.main()
