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
import espressomd.virtual_sites
import espressomd.integrate
from espressomd.shapes import Sphere, Wall

modes = {x for mode in set("@TEST_COMBINATION@".upper().split('-'))
         for x in [mode, mode.split('.')[0]]}

LB = ('LB.CPU' in modes or 'LB.GPU' in modes and espressomd.gpu_available())

EK = ('EK.GPU' in modes and espressomd.gpu_available()
      and espressomd.has_features('ELECTROKINETICS'))


class CheckpointTest(ut.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.checkpoint = espressomd.checkpointing.Checkpoint(
            checkpoint_id="mycheckpoint_@TEST_COMBINATION@_@TEST_BINARY@".replace(
                '.', '__'),
            checkpoint_path="@CMAKE_CURRENT_BINARY_DIR@")
        cls.checkpoint.load(0)

    @ut.skipIf(not LB, "Skipping test due to missing mode.")
    def test_LB(self):
        lbf = system.actors[0]
        cpt_mode = int("@TEST_BINARY@")
        cpt_path = self.checkpoint.checkpoint_dir + "/lb{}.cpt"
        with self.assertRaises(RuntimeError):
            lbf.load_checkpoint(cpt_path.format("-corrupted"), cpt_mode)
        with self.assertRaisesRegex(RuntimeError, 'grid dimensions mismatch'):
            lbf.load_checkpoint(cpt_path.format("-wrong-boxdim"), cpt_mode)
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
        reference = {'agrid': 0.5, 'visc': 1.3, 'dens': 1.5, 'tau': 0.01}
        for key in reference:
            self.assertIn(key, state)
            self.assertAlmostEqual(reference[key], state[key], delta=1E-7)

    @utx.skipIfMissingFeatures('ELECTROKINETICS')
    @ut.skipIf(not EK, "Skipping test due to missing mode.")
    def test_EK(self):
        ek = system.actors[0]
        ek_species = ek.get_params()['species'][0]
        cpt_path = self.checkpoint.checkpoint_dir + "/ek"
        ek.load_checkpoint(cpt_path)
        precision = 5
        m = np.pi / 12
        nx = int(np.round(system.box_l[0] / ek.get_params()["agrid"]))
        ny = int(np.round(system.box_l[1] / ek.get_params()["agrid"]))
        nz = int(np.round(system.box_l[2] / ek.get_params()["agrid"]))
        grid_3D = np.fromfunction(
            lambda i, j, k: np.cos(i * m) * np.cos(j * m) * np.cos(k * m),
            (nx, ny, nz), dtype=float)
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    np.testing.assert_almost_equal(
                        np.copy(ek_species[i, j, k].density),
                        grid_3D[i, j, k],
                        decimal=precision)
        state = ek.get_params()
        reference = {'agrid': 0.5, 'lb_density': 26.15,
                     'viscosity': 1.7, 'friction': 0.0,
                     'T': 1.1, 'prefactor': 0.88, 'stencil': "linkcentered"}
        for key in reference:
            self.assertIn(key, state)
            self.assertAlmostEqual(reference[key], state[key], delta=1E-5)
        state_species = ek_species.get_params()
        reference_species = {'density': 0.4, 'D': 0.02, 'valency': 0.3}
        for key in reference_species:
            self.assertIn(key, state_species)
            self.assertAlmostEqual(
                reference_species[key],
                state_species[key],
                delta=1E-5)
        self.assertAlmostEqual(
            state_species['ext_force_density'][0],
            0.01,
            delta=1E-5)
        self.assertAlmostEqual(
            state_species['ext_force_density'][1],
            -0.08,
            delta=1E-5)
        self.assertAlmostEqual(
            state_species['ext_force_density'][2],
            0.06,
            delta=1E-5)

    def test_variables(self):
        self.assertEqual(system.cell_system.skin, 0.1)
        self.assertEqual(system.time_step, 0.01)
        self.assertEqual(system.min_global_cut, 2.0)

    def test_part(self):
        np.testing.assert_allclose(
            np.copy(system.part[0].pos), np.array([1.0, 2.0, 3.0]))
        np.testing.assert_allclose(
            np.copy(system.part[1].pos), np.array([1.0, 1.0, 2.0]))
        np.testing.assert_allclose(np.copy(system.part[0].f), particle_force0)
        np.testing.assert_allclose(np.copy(system.part[1].f), particle_force1)

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
        self.assertEqual(thmst['seed'], 42 + 6)

    @utx.skipIfMissingFeatures('NPT')
    @ut.skipIf('THERM.NPT' not in modes, 'NPT thermostat not in modes')
    def test_thermostat_NPT(self):
        thmst = system.thermostat.get_state()[0]
        self.assertEqual(thmst['type'], 'NPT_ISO')
        self.assertEqual(thmst['seed'], 42)
        self.assertEqual(thmst['gamma0'], 2.0)
        self.assertEqual(thmst['gammav'], 0.1)

    @utx.skipIfMissingFeatures('NPT')
    @ut.skipIf('INT.NPT' not in modes, 'NPT integrator not in modes')
    def test_integrator_NPT(self):
        integ = system.integrator.get_state()
        self.assertIsInstance(
            integ, espressomd.integrate.VelocityVerletIsotropicNPT)
        params = integ.get_params()
        self.assertEqual(params['ext_pressure'], 2.0)
        self.assertEqual(params['piston'], 0.01)
        self.assertEqual(params['direction'], [1, 0, 0])
        self.assertEqual(params['cubic_box'], False)

    @ut.skipIf('INT.SD' not in modes, 'SD integrator not in modes')
    def test_integrator_SD(self):
        integ = system.integrator.get_state()
        self.assertIsInstance(integ, espressomd.integrate.SteepestDescent)
        params = integ.get_params()
        self.assertEqual(params['f_max'], 2.0)
        self.assertEqual(params['gamma'], 0.1)
        self.assertEqual(params['max_displacement'], 0.01)

    @ut.skipIf('INT.NVT' not in modes, 'NVT integrator not in modes')
    def test_integrator_NVT(self):
        integ = system.integrator.get_state()
        self.assertIsInstance(integ, espressomd.integrate.VelocityVerlet)
        params = integ.get_params()
        self.assertEqual(params, {})

    @ut.skipIf('INT' in modes, 'VV integrator not the default')
    def test_integrator_VV(self):
        integ = system.integrator.get_state()
        self.assertIsInstance(integ, espressomd.integrate.VelocityVerlet)
        params = integ.get_params()
        self.assertEqual(params, {})

    @ut.skipIf('INT.BD' not in modes, 'BD integrator not in modes')
    def test_integrator_BD(self):
        integ = system.integrator.get_state()
        self.assertIsInstance(integ, espressomd.integrate.BrownianDynamics)
        params = integ.get_params()
        self.assertEqual(params, {})

    @utx.skipIfMissingFeatures('LENNARD_JONES')
    @ut.skipIf('LJ' not in modes, "Skipping test due to missing mode.")
    def test_non_bonded_inter(self):
        state = system.non_bonded_inter[
            0, 0].lennard_jones._get_params_from_es_core()
        state2 = system.non_bonded_inter[
            3, 0].lennard_jones._get_params_from_es_core()
        reference = {'shift': 0.1, 'sigma': 1.3, 'epsilon': 1.2,
                     'cutoff': 2.0, 'offset': 0.0, 'min': 0.0}
        reference2 = {'shift': 0.1, 'sigma': 1.7, 'epsilon': 1.2,
                      'cutoff': 2.0, 'offset': 0.0, 'min': 0.0}
        self.assertEqual(
            len(set(state.items()) & set(reference.items())), len(reference))
        self.assertEqual(len(set(state2.items()) & set(
            reference2.items())), len(reference2))

    def test_bonded_inter(self):
        state = system.part[1].bonds[0][0].params
        reference = {'r_0': 0.0, 'k': 1.0}
        for key in reference:
            self.assertAlmostEqual(state[key], reference[key], delta=1E-10)
        if 'THERM.LB' not in modes:
            state = system.part[1].bonds[1][0].params
            reference = {'temp_com': 0., 'gamma_com': 0., 'temp_distance': 0.2,
                         'gamma_distance': 0.5, 'r_cut': 2.0, 'seed': 51}
            for key in reference:
                self.assertAlmostEqual(state[key], reference[key], delta=1E-10)

    @utx.skipIfMissingFeatures(['VIRTUAL_SITES', 'VIRTUAL_SITES_RELATIVE'])
    def test_virtual_sites(self):
        self.assertTrue(system.part[1].virtual)
        self.assertIsInstance(
            system.virtual_sites,
            espressomd.virtual_sites.VirtualSitesRelative)

    def test_mean_variance_calculator(self):
        np.testing.assert_array_equal(
            acc.get_mean(), np.array([[1.0, 1.5, 2.0], [1.0, 1.0, 2.0]]))
        np.testing.assert_array_equal(
            acc.get_variance(), np.array([[0., 0.5, 2.], [0., 0., 0.]]))
        np.testing.assert_array_equal(
            system.auto_update_accumulators[0].get_variance(),
            np.array([[0., 0.5, 2.], [0., 0., 0.]]))

    @utx.skipIfMissingFeatures('P3M')
    @ut.skipIf('P3M.CPU' not in modes,
               "Skipping test due to missing combination.")
    def test_p3m(self):
        self.assertTrue(any(isinstance(actor, espressomd.electrostatics.P3M)
                            for actor in system.actors.active_actors))

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

    @ut.skipIf(not LB or EK or not (espressomd.has_features("LB_BOUNDARIES")
                                    or espressomd.has_features("LB_BOUNDARIES_GPU")), "Missing features")
    def test_lb_boundaries(self):
        self.assertEqual(len(system.lbboundaries), 1)
        np.testing.assert_allclose(
            np.copy(system.lbboundaries[0].velocity), [1e-4, 1e-4, 0])
        self.assertEqual(type(system.lbboundaries[0].shape), Wall)

    def test_constraints(self):
        from espressomd import constraints
        self.assertEqual(len(system.constraints),
                         8 - int(not espressomd.has_features("ELECTROSTATICS")))
        c = system.constraints

        self.assertEqual(type(c[0].shape), Sphere)
        self.assertAlmostEqual(c[0].shape.radius, 0.1, delta=1E-10)
        self.assertEqual(c[0].particle_type, 17)

        self.assertEqual(type(c[1].shape), Wall)
        np.testing.assert_allclose(np.copy(c[1].shape.normal),
                                   [1. / np.sqrt(3)] * 3)

        self.assertEqual(type(c[2]), constraints.Gravity)
        np.testing.assert_allclose(np.copy(c[2].g), [1., 2., 3.])

        self.assertEqual(type(c[3]), constraints.HomogeneousMagneticField)
        np.testing.assert_allclose(np.copy(c[3].H), [1., 2., 3.])

        self.assertEqual(type(c[4]), constraints.HomogeneousFlowField)
        np.testing.assert_allclose(np.copy(c[4].u), [1., 2., 3.])
        self.assertAlmostEqual(c[4].gamma, 2.3, delta=1E-10)

        self.assertEqual(type(c[5]), constraints.PotentialField)
        self.assertEqual(c[5].field.shape, (14, 16, 18, 1))
        self.assertAlmostEqual(c[5].default_scale, 1.6, delta=1E-10)
        np.testing.assert_allclose(np.copy(c[5].origin), [-0.5, -0.5, -0.5])
        np.testing.assert_allclose(np.copy(c[5].grid_spacing), np.ones(3))
        ref_pot = constraints.PotentialField(
            field=pot_field_data, grid_spacing=np.ones(3), default_scale=1.6)
        np.testing.assert_allclose(np.copy(c[5].field), np.copy(ref_pot.field),
                                   atol=1e-10)

        self.assertEqual(type(c[6]), constraints.ForceField)
        self.assertEqual(c[6].field.shape, (14, 16, 18, 3))
        self.assertAlmostEqual(c[6].default_scale, 1.4, delta=1E-10)
        np.testing.assert_allclose(np.copy(c[6].origin), [-0.5, -0.5, -0.5])
        np.testing.assert_allclose(np.copy(c[6].grid_spacing), np.ones(3))
        ref_vec = constraints.ForceField(
            field=vec_field_data, grid_spacing=np.ones(3), default_scale=1.4)
        np.testing.assert_allclose(np.copy(c[6].field), np.copy(ref_vec.field),
                                   atol=1e-10)

        if espressomd.has_features("ELECTROSTATICS"):
            self.assertEqual(type(c[7]), constraints.ElectricPlaneWave)
            np.testing.assert_allclose(np.copy(c[7].E0), [1., -2., 3.])
            np.testing.assert_allclose(np.copy(c[7].k), [-.1, .2, .3])
            self.assertAlmostEqual(c[7].omega, 5., delta=1E-10)
            self.assertAlmostEqual(c[7].phi, 1.4, delta=1E-10)


if __name__ == '__main__':
    ut.main()
