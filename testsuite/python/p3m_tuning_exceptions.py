#
# Copyright (C) 2020-2022 The ESPResSo project
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
import espressomd
import espressomd.utils
import espressomd.electrostatics
import espressomd.magnetostatics
import unittest as ut
import unittest_decorators as utx
import itertools
import numpy as np


class Test(ut.TestCase):

    system = espressomd.System(box_l=[10., 10., 10.])

    # tune parameters that are valid for all CPU methods
    valid_params = {
        'P3MGPU':
        {'cao': 2, 'r_cut': 3.1836, 'accuracy': 0.01, 'mesh': [8, 8, 8],
         'mesh_off': [0.5, 0.5, 0.5], 'prefactor': 2.0, 'alpha': 0.5115},
        'P3M':
        {'cao': 2, 'r_cut': 3.1836, 'accuracy': 0.01, 'mesh': [8, 8, 8],
         'mesh_off': [0.5, 0.5, 0.5], 'prefactor': 2.0, 'alpha': 0.5115},
        'DP3M':
        {'cao': 1, 'r_cut': 3.28125, 'accuracy': 0.01, 'mesh': [5, 5, 5],
         'mesh_off': [0.5, 0.5, 0.5], 'prefactor': 2.0, 'alpha': 0.455},
    }

    def get_valid_params(self, key, **custom_params):
        params = self.valid_params[key].copy()
        params.update(custom_params)
        return params

    def setUp(self):
        self.system.box_l = [10., 10., 10.]
        self.system.periodicity = [True, True, True]

    def tearDown(self):
        self.system.actors.clear()
        self.system.part.clear()
        # assertion: there should be no pending runtime error
        espressomd.utils.handle_errors("tearDown")

    def add_charged_particles(self):
        self.system.part.add(pos=[[0., 0., 0.], [0.5, 0.5, 0.5]], q=[-1., 1.])

    def add_magnetic_particles(self):
        self.system.part.add(pos=[[0.01, 0.01, 0.01], [0.5, 0.5, 0.5]],
                             dip=[(1., 0., 0.), (-1., 0., 0.)],
                             rotation=2 * [(True, True, True)])

    ##################################################
    # block of tests where the time_step is negative #
    ##################################################

    @utx.skipIfMissingGPU()
    @utx.skipIfMissingFeatures("P3M")
    def test_01_time_not_set_p3m_gpu(self):
        self.add_charged_particles()

        solver = espressomd.electrostatics.P3MGPU(prefactor=2, accuracy=1e-2)
        with self.assertRaisesRegex(Exception, 'time_step not set'):
            self.system.actors.add(solver)

    @utx.skipIfMissingFeatures("P3M")
    def test_01_time_not_set_p3m_cpu(self):
        self.add_charged_particles()

        solver = espressomd.electrostatics.P3M(prefactor=2, accuracy=1e-2)
        with self.assertRaisesRegex(Exception, 'time_step not set'):
            self.system.actors.add(solver)

    @utx.skipIfMissingFeatures("DP3M")
    def test_01_time_not_set_dp3m_cpu(self):
        self.add_magnetic_particles()

        solver = espressomd.magnetostatics.DipolarP3M(
            prefactor=2, accuracy=1e-2)
        with self.assertRaisesRegex(Exception, 'time_step not set'):
            self.system.actors.add(solver)

    ##############################################
    # block of tests where particles are missing #
    ##############################################

    @utx.skipIfMissingGPU()
    @utx.skipIfMissingFeatures("P3M")
    def test_02_no_particles_p3m_gpu(self):
        self.system.time_step = 0.01

        solver = espressomd.electrostatics.P3MGPU(prefactor=2, accuracy=1e-2)
        with self.assertRaisesRegex(RuntimeError, 'no charged particles in the system'):
            self.system.actors.add(solver)

    @utx.skipIfMissingFeatures("P3M")
    def test_02_no_particles_p3m_cpu(self):
        self.system.time_step = 0.01

        solver = espressomd.electrostatics.P3M(prefactor=2, accuracy=1e-2)
        with self.assertRaisesRegex(RuntimeError, 'no charged particles in the system'):
            self.system.actors.add(solver)

    @utx.skipIfMissingFeatures("DP3M")
    def test_02_no_particles_dp3m_cpu(self):
        self.system.time_step = 0.01

        solver = espressomd.magnetostatics.DipolarP3M(
            **self.get_valid_params('DP3M'))
        with self.assertRaisesRegex(RuntimeError, 'DipolarP3M: no dipolar particles in the system'):
            self.system.actors.add(solver)

    @utx.skipIfMissingFeatures("DP3M")
    def test_02_accuracy_dp3m_cpu(self):
        self.system.time_step = 0.01
        self.add_magnetic_particles()

        solver = espressomd.magnetostatics.DipolarP3M(
            **self.get_valid_params('DP3M', accuracy=1e-20))
        with self.assertRaisesRegex(Exception, 'DipolarP3M: failed to reach requested accuracy'):
            self.system.actors.add(solver)

    #######################################
    # block of tests with non-cubic boxes #
    #######################################

    @utx.skipIfMissingGPU()
    @utx.skipIfMissingFeatures("P3M")
    def test_03_non_cubic_box_p3m_gpu(self):
        self.system.box_l = [10., 10., 20.]
        self.system.time_step = 0.01
        self.add_charged_particles()

        solver = espressomd.electrostatics.P3MGPU(
            prefactor=2, accuracy=1e-2, epsilon=1)
        with self.assertRaisesRegex(RuntimeError, 'P3M: non-metallic epsilon requires cubic box'):
            self.system.actors.add(solver)

    @utx.skipIfMissingFeatures("P3M")
    def test_03_non_cubic_box_p3m_cpu(self):
        self.system.box_l = [10., 10., 20.]
        self.system.time_step = 0.01
        self.add_charged_particles()

        solver = espressomd.electrostatics.P3M(
            prefactor=2, accuracy=1e-2, epsilon=1, mesh=[8, 8, 8])
        with self.assertRaisesRegex(RuntimeError, 'P3M: non-metallic epsilon requires cubic box'):
            self.system.actors.add(solver)

        self.system.box_l = [10., 10., 10.]
        solver = espressomd.electrostatics.P3M(
            prefactor=2, accuracy=1e-2, epsilon=1, mesh=[4, 8, 8])
        with self.assertRaisesRegex(RuntimeError, 'P3M: non-metallic epsilon requires cubic box'):
            self.system.actors.add(solver)

    @utx.skipIfMissingFeatures("DP3M")
    def test_03_non_cubic_box_dp3m_cpu(self):
        self.system.box_l = [10., 10., 20.]
        self.system.time_step = 0.01
        self.add_magnetic_particles()

        with self.assertRaisesRegex(RuntimeError, 'DipolarP3M: requires a cubic box'):
            self.system.actors.add(espressomd.magnetostatics.DipolarP3M(
                **self.get_valid_params('DP3M'), tune=False))

    ##########################################
    # block of tests with invalid parameters #
    ##########################################

    def check_invalid_params(self, class_solver, **custom_params):
        valid_params = {
            'prefactor': 2.0, 'accuracy': .01, 'tune': False, 'cao': 3,
            'r_cut': 0.373, 'alpha': 3.81, 'mesh': (8, 8, 8), 'epsilon': 0.,
            'mesh_off': [0.5, 0.5, 0.5]}
        valid_params.update(custom_params)

        invalid_params = [
            ('prefactor', -2.0, "Parameter 'prefactor' must be > 0"),
            ('epsilon', -1.0, "Parameter 'epsilon' must be >= 0"),
            ('cao', 0, "Parameter 'cao' must be >= 1 and <= 7"),
            ('cao', 8, "Parameter 'cao' must be >= 1 and <= 7"),
            ('r_cut', -2.0, "Parameter 'r_cut' must be > 0"),
            ('alpha', -2.0, "Parameter 'alpha' must be > 0"),
            ('accuracy', -2.0, "Parameter 'accuracy' must be > 0"),
            ('mesh', (-1, -1, -1), "Parameter 'mesh' must be > 0"),
            ('mesh', (2, 2, 2), "Parameter 'cao' cannot be larger than 'mesh'"),
            ('mesh_off', (-2, 1, 1), "Parameter 'mesh_off' must be >= 0 and <= 1"),
        ]
        if class_solver is espressomd.magnetostatics.DipolarP3M:
            invalid_params.append(
                ('mesh', (4, 2, 2), "DipolarP3M requires a cubic mesh")
            )

        for key, invalid_value, err_msg in invalid_params:
            params = valid_params.copy()
            params[key] = invalid_value
            with self.assertRaisesRegex(ValueError, err_msg):
                class_solver(**params)

        # cannot add an actor if cell system isn't compatible
        with self.assertRaisesRegex(RuntimeError, "P3M: requires periodicity"):
            self.system.periodicity = (True, True, False)
            solver = class_solver(**valid_params)
            self.system.actors.add(solver)
        self.assertEqual(len(self.system.actors), 0)
        self.system.periodicity = (True, True, True)

    def check_invalid_cell_systems(self):
        # check periodicity exceptions
        for periodicity in itertools.product(range(2), range(2), range(2)):
            if periodicity == (1, 1, 1):
                continue
            with self.assertRaisesRegex(Exception, r"P3M: requires periodicity \(1 1 1\)"):
                self.system.periodicity = periodicity
        self.system.periodicity = (1, 1, 1)

        # check cell system exceptions
        with self.assertRaisesRegex(Exception, "P3M: requires the regular or hybrid decomposition cell system"):
            self.system.cell_system.set_n_square()
            self.system.analysis.energy()
        self.system.cell_system.set_regular_decomposition()

    @utx.skipIfMissingFeatures("P3M")
    def test_04_invalid_params_p3m_cpu(self):
        self.system.time_step = 0.01
        self.add_charged_particles()

        self.check_invalid_params(espressomd.electrostatics.P3M)

        # set up a valid actor
        solver = espressomd.electrostatics.P3M(
            prefactor=2, accuracy=0.1, cao=2, r_cut=3.18, mesh=8)
        self.system.actors.add(solver)
        self.check_invalid_cell_systems()

    @utx.skipIfMissingGPU()
    @utx.skipIfMissingFeatures("P3M")
    def test_04_invalid_params_p3m_gpu(self):
        self.system.time_step = 0.01
        self.add_charged_particles()

        self.check_invalid_params(espressomd.electrostatics.P3MGPU,
                                  mesh=3 * [28], alpha=0.3548, r_cut=4.4434)

    @utx.skipIfMissingFeatures("DP3M")
    def test_04_invalid_params_dp3m_cpu(self):
        self.system.time_step = 0.01
        self.add_magnetic_particles()

        self.check_invalid_params(espressomd.magnetostatics.DipolarP3M)

        # check bisection exception
        with self.assertRaisesRegex(RuntimeError, r"Root must be bracketed for bisection in dp3m_rtbisection"):
            solver = espressomd.magnetostatics.DipolarP3M(
                prefactor=2, accuracy=0.01, cao=1,
                r_cut=0.373, alpha=3.81, mesh=(8, 8, 8))
            self.system.actors.add(solver)
        self.system.actors.clear()

        # set up a valid actor
        solver = espressomd.magnetostatics.DipolarP3M(
            **self.get_valid_params('DP3M'), tune=False)
        self.system.actors.add(solver)
        self.check_invalid_cell_systems()

    def check_invalid_params_layer_corrections(self, solver_p3m, class_lc):
        with self.assertRaisesRegex(ValueError, "Parameter 'gap_size' must be > 0"):
            class_lc(actor=solver_p3m, gap_size=-1., maxPWerror=0.01)
        with self.assertRaisesRegex(ValueError, "Parameter 'maxPWerror' must be > 0"):
            class_lc(actor=solver_p3m, gap_size=1., maxPWerror=0.)
        with self.assertRaisesRegex(ValueError, "Parameter 'far_cut' must be > 0"):
            class_lc(actor=solver_p3m, gap_size=1., maxPWerror=1., far_cut=0.)
        with self.assertRaisesRegex(ValueError, "Parameter 'far_cut' must be > 0"):
            class_lc(actor=solver_p3m, gap_size=1., maxPWerror=1., far_cut=-2.)
        with self.assertRaisesRegex(RuntimeError, "LC gap size .+ larger than box length in z-direction"):
            lc = class_lc(actor=solver_p3m, gap_size=100., maxPWerror=0.01)
            self.system.actors.add(lc)
        self.assertEqual(len(self.system.actors), 0)

    def check_invalid_params_elc_p3m(self, solver_p3m):
        ELC = espressomd.electrostatics.ELC
        self.check_invalid_params_layer_corrections(solver_p3m, ELC)

        self.system.part.by_id(0).q = -1.00001

        with self.assertRaisesRegex(RuntimeError, "ELC does not currently support non-neutral systems with a dielectric contrast"):
            actor = ELC(actor=solver_p3m, gap_size=1., maxPWerror=1.,
                        pot_diff=3., delta_mid_top=0.5, delta_mid_bot=0.5,
                        const_pot=True, check_neutrality=False)
            self.system.actors.add(actor)
        self.system.actors.clear()

        with self.assertRaisesRegex(RuntimeError, "ELC does not work for non-neutral systems and non-metallic dielectric contrast"):
            actor = ELC(actor=solver_p3m, gap_size=1., maxPWerror=1.,
                        pot_diff=0., delta_mid_top=0.5, delta_mid_bot=0.5,
                        const_pot=False, check_neutrality=False)
            self.system.actors.add(actor)
        self.system.actors.clear()

        self.system.part.by_id(0).q = -1

        with self.assertRaisesRegex(RuntimeError, "ELC tuning failed: maxPWerror too small"):
            # reduce box size to make tuning converge in at most 50 steps
            self.system.box_l = [1., 1., 1.]
            elc = ELC(actor=solver_p3m, gap_size=0.5, maxPWerror=1e-90)
            self.system.actors.add(elc)
        self.assertEqual(len(self.system.actors), 0)
        self.system.box_l = [10., 10., 10.]

        # r_cut > gap isn't allowed with dielectric contrasts
        p3m = espressomd.electrostatics.P3M(
            prefactor=1.5, r_cut=3., accuracy=0.01, mesh=[8, 8, 8])
        with self.assertRaisesRegex(RuntimeError, "failed to reach requested accuracy"):
            elc = ELC(actor=p3m, gap_size=p3m.r_cut / 2., maxPWerror=0.01,
                      delta_mid_top=0.5, delta_mid_bot=0.5, pot_diff=-3.,
                      const_pot=True)
            self.system.actors.add(elc)
        self.assertEqual(len(self.system.actors), 0)

        # r_cut > gap is allowed without dielectric contrasts
        elc = ELC(actor=p3m, gap_size=p3m.r_cut / 2., maxPWerror=0.01)
        self.system.actors.add(elc)
        self.assertEqual(len(self.system.actors), 1)
        self.assertAlmostEqual(elc.prefactor, 1.5, delta=1e-12)

    @utx.skipIfMissingFeatures("P3M")
    def test_04_invalid_params_elc_p3m_cpu(self):
        self.system.time_step = 0.01
        self.add_charged_particles()
        params = self.get_valid_params('P3M', tune=False)
        solver = espressomd.electrostatics.P3M(**params)
        self.check_invalid_params_elc_p3m(solver)

    @utx.skipIfMissingGPU()
    @utx.skipIfMissingFeatures("P3M")
    def test_04_invalid_params_elc_p3m_gpu(self):
        self.system.time_step = 0.01
        self.add_charged_particles()
        params = self.get_valid_params('P3MGPU', tune=False)
        solver = espressomd.electrostatics.P3MGPU(**params)
        self.check_invalid_params_elc_p3m(solver)

    @utx.skipIfMissingFeatures("DP3M")
    def test_04_invalid_params_dlc_dp3m_cpu(self):
        self.system.time_step = 0.01
        self.add_magnetic_particles()

        dp3m_params = {'accuracy': 1e-6, 'mesh': [25, 25, 25], 'cao': 7,
                       'prefactor': 1.1, 'r_cut': 4.50, 'alpha': 0.8216263}

        solver_dp3m = espressomd.magnetostatics.DipolarP3M(
            epsilon='metallic', tune=False, **dp3m_params)
        self.check_invalid_params_layer_corrections(
            solver_dp3m, espressomd.magnetostatics.DLC)

        solver_mdlc = espressomd.magnetostatics.DLC(
            gap_size=1, maxPWerror=1e-30, actor=solver_dp3m)
        with self.assertRaisesRegex(RuntimeError, "DLC tuning failed: maxPWerror too small"):
            self.system.actors.add(solver_mdlc)

        dp3m_params = {'accuracy': 1e-30, 'mesh': [6, 6, 6],
                       'prefactor': 1.1, 'r_cut': 4.50}

        solver_dp3m = espressomd.magnetostatics.DipolarP3M(
            prefactor=1., accuracy=1e-30, mesh=[6, 6, 6], r_cut=4.50)
        solver_mdlc = espressomd.magnetostatics.DLC(
            gap_size=1., maxPWerror=1e-2, actor=solver_dp3m)
        with self.assertRaisesRegex(RuntimeError, "P3M: failed to reach requested accuracy"):
            self.system.actors.add(solver_mdlc)

    ###########################################################
    # block of tests where tuning should not throw exceptions #
    ###########################################################

    @utx.skipIfMissingGPU()
    @utx.skipIfMissingFeatures("P3M")
    def test_09_no_errors_p3m_gpu(self):
        self.system.time_step = 0.01
        self.add_charged_particles()

        # mesh is fixed to significantly speed up tuning
        solver = espressomd.electrostatics.P3MGPU(
            prefactor=2, accuracy=1e-2, epsilon='metallic', mesh=[20, 20, 20])
        self.system.actors.add(solver)

    @utx.skipIfMissingFeatures("P3M")
    def test_09_no_errors_p3m_cpu(self):
        self.system.time_step = 0.01
        self.add_charged_particles()

        valid_params = {
            'mesh_off': [-1., -1., -1.],  # sentinel
            'cao': 2, 'r_cut': 3.18, 'mesh': 8}

        # tuning with cao or r_cut or mesh constrained, or without constraints
        for key, value in valid_params.items():
            solver = espressomd.electrostatics.P3M(
                prefactor=2, accuracy=1e-2, epsilon=0.0, **{key: value})
            self.system.actors.add(solver)
            self.system.actors.clear()

    @utx.skipIfMissingFeatures("DP3M")
    def test_09_no_errors_dp3m_cpu(self):
        self.system.time_step = 0.01
        self.add_magnetic_particles()

        valid_params = {
            'mesh_off': [-1., -1., -1.],  # sentinel
            'cao': 1, 'r_cut': 3.28125, 'mesh': 5}

        # tuning with cao or r_cut or mesh constrained, or without constraints
        for key, value in valid_params.items():
            solver = espressomd.magnetostatics.DipolarP3M(
                prefactor=2, accuracy=1e-2, **{key: value})
            self.system.actors.add(solver)
            self.system.actors.clear()

    @utx.skipIfMissingFeatures("P3M")
    def test_09_no_errors_p3m_cpu_rescale_mesh(self):
        self.system.box_l = [10., 15., 20.]
        self.system.time_step = 0.01
        self.add_charged_particles()

        solver = espressomd.electrostatics.P3M(prefactor=2, accuracy=1e-2,
                                               epsilon='metallic',
                                               mesh=[8, -1, -1])
        self.system.actors.add(solver)
        self.assertEqual(solver.mesh, [8, 12, 16])

        # check MD cell reset event
        self.system.box_l = self.system.box_l
        self.system.periodicity = self.system.periodicity
        self.system.cell_system.node_grid = self.system.cell_system.node_grid

    @utx.skipIfMissingGPU()
    @utx.skipIfMissingFeatures("P3M")
    def test_09_no_errors_p3m_gpu_rescale_mesh(self):
        self.system.box_l = [10., 10., 20.]
        self.system.time_step = 0.01
        self.add_charged_particles()

        solver = espressomd.electrostatics.P3MGPU(prefactor=2, accuracy=1e-1,
                                                  epsilon='metallic',
                                                  mesh=[20, -1, -1])
        self.system.actors.add(solver)
        self.assertEqual(solver.mesh, [20, 20, 40])

        # check MD cell reset event
        self.system.box_l = self.system.box_l
        self.system.periodicity = self.system.periodicity
        self.system.cell_system.node_grid = self.system.cell_system.node_grid

    @utx.skipIfMissingFeatures("DP3M")
    def test_09_no_errors_dp3m_cpu_rescale_mesh(self):
        self.system.time_step = 0.01
        self.add_magnetic_particles()

        dp3m_params = {'accuracy': 1e-6, 'mesh': [25, 25, 25], 'cao': 7,
                       'prefactor': 1.1, 'r_cut': 4.50, 'alpha': 0.8216263}
        solver = espressomd.magnetostatics.DipolarP3M(
            epsilon='metallic', tune=False, **dp3m_params)
        self.system.actors.add(solver)

        # check MD cell reset event
        self.system.box_l = self.system.box_l
        self.system.periodicity = self.system.periodicity
        self.system.cell_system.node_grid = self.system.cell_system.node_grid

    def check_tuning_layer_corrections(self, class_p3m, class_lc, params):
        if class_p3m is espressomd.magnetostatics.DipolarP3M:
            mesh_a = np.array([2., 2., 2.])
        else:
            mesh_a = np.array([2., 4., 8.])
        self.system.box_l = mesh_a * params["mesh"]
        self.system.time_step = 0.01
        non_metallic_epsilon = 20.
        p3m = class_p3m(epsilon=non_metallic_epsilon, **params)
        self.assertEqual(p3m.epsilon, non_metallic_epsilon)
        lc = class_lc(actor=p3m, gap_size=1., maxPWerror=1.)
        # Non-metallic epsilon values are not allowed for non-cubic boxes and
        # will cause tuning to fail. Since ELC doesn't support non-metallic
        # epsilon values either, it sets metallic epsilon before tuning.
        if class_lc is espressomd.electrostatics.ELC:
            self.assertEqual(p3m.epsilon, 0.)
        self.system.actors.add(lc)

        # check parameter rescaling
        alpha = p3m.alpha
        r_cut = params["r_cut"]
        r_cut_iL = r_cut / self.system.box_l[0]
        alpha_L = alpha * self.system.box_l[0]
        np.testing.assert_allclose(np.copy(p3m.a), mesh_a, atol=1e-12)
        np.testing.assert_allclose(p3m.r_cut, r_cut, atol=1e-12)
        np.testing.assert_allclose(p3m.r_cut_iL, r_cut_iL, atol=1e-12)
        np.testing.assert_allclose(p3m.alpha_L, alpha_L, atol=1e-12)
        mesh_a = np.array([4., 4., 4.])
        self.system.box_l = mesh_a * params["mesh"]
        np.testing.assert_allclose(np.copy(p3m.a), mesh_a, atol=1e-12)
        np.testing.assert_allclose(p3m.r_cut, r_cut * 2., atol=1e-12)
        np.testing.assert_allclose(p3m.r_cut_iL, r_cut_iL, atol=1e-12)
        np.testing.assert_allclose(p3m.alpha, alpha / 2., atol=1e-12)
        np.testing.assert_allclose(p3m.alpha_L, alpha_L, atol=1e-12)

    @utx.skipIfMissingFeatures("P3M")
    def test_09_no_errors_elc_p3m_cpu_rescale_mesh(self):
        self.add_charged_particles()
        self.check_tuning_layer_corrections(
            espressomd.electrostatics.P3M,
            espressomd.electrostatics.ELC,
            self.get_valid_params("P3M", accuracy=0.1))

    @utx.skipIfMissingGPU()
    @utx.skipIfMissingFeatures("P3M")
    def test_09_no_errors_elc_p3m_gpu_rescale_mesh(self):
        self.add_charged_particles()
        self.check_tuning_layer_corrections(
            espressomd.electrostatics.P3MGPU,
            espressomd.electrostatics.ELC,
            self.get_valid_params("P3MGPU", accuracy=0.1))

    @utx.skipIfMissingFeatures("DP3M")
    def test_09_no_errors_dlc_dp3m_cpu_rescale_mesh(self):
        self.add_magnetic_particles()
        self.check_tuning_layer_corrections(
            espressomd.magnetostatics.DipolarP3M,
            espressomd.magnetostatics.DLC,
            self.get_valid_params("DP3M", accuracy=0.1))


if __name__ == "__main__":
    ut.main(failfast=True)
