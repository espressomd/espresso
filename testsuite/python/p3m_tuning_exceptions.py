#
# Copyright (C) 2020 The ESPResSo project
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
import unittest as ut
import unittest_decorators as utx


class P3M_tuning_test(ut.TestCase):

    system = espressomd.System(box_l=[10., 10., 10.])

    def setUp(self):
        self.system.box_l = [10., 10., 10.]
        self.system.periodicity = [True, True, True]

    def tearDown(self):
        self.system.actors.clear()
        self.system.part.clear()

    def add_charged_particles(self):
        self.system.part.add(pos=[[0, 0, 0], [.5, .5, .5]], q=[-1, 1])

    def add_magnetic_particles(self):
        self.system.part.add(pos=[[0.01, 0.01, 0.01], [.5, .5, .5]],
                             rotation=2 * [(1, 1, 1)], dip=2 * [(1, 0, 0)])

    ##################################################
    # block of tests where the time_step is negative #
    ##################################################

    @utx.skipIfMissingGPU()
    @utx.skipIfMissingFeatures("P3M")
    def test_01_time_not_set_p3m_gpu(self):
        import espressomd.electrostatics

        self.add_charged_particles()

        solver = espressomd.electrostatics.P3MGPU(prefactor=2, accuracy=1e-2)
        with self.assertRaisesRegex(Exception, 'P3M: tuning failed: ERROR: time_step not set'):
            self.system.actors.add(solver)

    @utx.skipIfMissingFeatures("P3M")
    def test_01_time_not_set_p3m_cpu(self):
        import espressomd.electrostatics

        self.add_charged_particles()

        solver = espressomd.electrostatics.P3M(prefactor=2, accuracy=1e-2)
        with self.assertRaisesRegex(Exception, 'P3M: tuning failed: ERROR: time_step not set'):
            self.system.actors.add(solver)

    @utx.skipIfMissingFeatures("DP3M")
    def test_01_time_not_set_dp3m_cpu(self):
        import espressomd.magnetostatics

        self.add_magnetic_particles()

        solver = espressomd.magnetostatics.DipolarP3M(
            prefactor=2, accuracy=1e-2)
        with self.assertRaisesRegex(Exception, 'P3M: tuning failed: ERROR: time_step not set'):
            self.system.actors.add(solver)

    ##############################################
    # block of tests where particles are missing #
    ##############################################

    @utx.skipIfMissingGPU()
    @utx.skipIfMissingFeatures("P3M")
    def test_02_no_particles_p3m_gpu(self):
        import espressomd.electrostatics

        self.system.time_step = 0.01

        solver = espressomd.electrostatics.P3MGPU(prefactor=2, accuracy=1e-2)
        with self.assertRaisesRegex(Exception, 'P3M: tuning failed: ERROR: no charged particles in the system'):
            self.system.actors.add(solver)

    @utx.skipIfMissingFeatures("P3M")
    def test_02_no_particles_p3m_cpu(self):
        import espressomd.electrostatics

        self.system.time_step = 0.01

        solver = espressomd.electrostatics.P3M(prefactor=2, accuracy=1e-2)
        with self.assertRaisesRegex(Exception, 'P3M: tuning failed: ERROR: no charged particles in the system'):
            self.system.actors.add(solver)

    @utx.skipIfMissingFeatures("DP3M")
    def test_02_no_particles_dp3m_cpu(self):
        import espressomd.magnetostatics

        self.system.time_step = 0.01

        solver = espressomd.magnetostatics.DipolarP3M(
            prefactor=2, accuracy=1e-2)
        with self.assertRaisesRegex(Exception, 'P3M: tuning failed: ERROR: no dipolar particles in the system'):
            self.system.actors.add(solver)

    #######################################
    # block of tests with non-cubic boxes #
    #######################################

    @utx.skipIfMissingGPU()
    @utx.skipIfMissingFeatures("P3M")
    def test_03_non_cubic_box_p3m_gpu(self):
        import espressomd.electrostatics

        self.system.box_l = [10., 10., 20.]
        self.system.time_step = 0.01
        self.add_charged_particles()

        solver = espressomd.electrostatics.P3MGPU(
            prefactor=2, accuracy=1e-2, epsilon=1)
        with self.assertRaisesRegex(Exception, 'P3M: tuning failed: ERROR: non-metallic epsilon requires cubic box'):
            self.system.actors.add(solver)

    @utx.skipIfMissingFeatures("P3M")
    def test_03_non_cubic_box_p3m_cpu(self):
        import espressomd.electrostatics

        self.system.box_l = [10., 10., 20.]
        self.system.time_step = 0.01
        self.add_charged_particles()

        solver = espressomd.electrostatics.P3M(
            prefactor=2, accuracy=1e-2, epsilon=1)
        with self.assertRaisesRegex(Exception, 'P3M: tuning failed: ERROR: non-metallic epsilon requires cubic box'):
            self.system.actors.add(solver)

    @utx.skipIfMissingFeatures("DP3M")
    def test_03_non_cubic_box_dp3m_cpu(self):
        import espressomd.magnetostatics

        self.system.box_l = [10., 10., 20.]
        self.system.time_step = 0.01
        self.add_magnetic_particles()

        solver = espressomd.magnetostatics.DipolarP3M(
            prefactor=2, accuracy=1e-2)
        with self.assertRaisesRegex(Exception, 'P3M: tuning failed: ERROR: dipolar P3M requires a cubic box'):
            self.system.actors.add(solver)

    ##########################################
    # block of tests with invalid parameters #
    ##########################################

    def check_invalid_params(self, solver_class, **custom_params):
        valid_params = {
            'prefactor': 2, 'accuracy': .01, 'tune': False, 'cao': 1,
            'r_cut': 0.373, 'alpha': 3.81, 'mesh': (8, 8, 8),
            'mesh_off': [-1, -1, -1]}
        valid_params.update(custom_params)

        invalid_params = [
            ('cao', 0, 'P3M: invalid cao'),
            ('cao', 8, 'P3M: invalid cao'),
            ('r_cut', -2.0, 'P3M: invalid r_cut'),
            ('alpha', -2.0, 'P3M: invalid alpha'),
            ('accuracy', -2.0, 'P3M: invalid accuracy'),
            ('mesh', (-1, -1, -1), 'P3M: invalid mesh size'),
            ('mesh', (0, 0, 0), 'P3M: cao larger than mesh size'),
            ('mesh_off', (-2, 1, 1), 'P3M: invalid mesh offset'),
        ]

        for key, invalid_value, err_msg in invalid_params:
            params = valid_params.copy()
            params[key] = invalid_value
            solver = solver_class(**params)
            with self.assertRaisesRegex(RuntimeError, err_msg):
                self.system.actors.add(solver)
            self.system.actors.clear()

    @utx.skipIfMissingFeatures("P3M")
    def test_04_invalid_params_p3m_cpu(self):
        import espressomd.electrostatics

        self.system.time_step = 0.01
        self.add_charged_particles()

        self.check_invalid_params(espressomd.electrostatics.P3M)

    @utx.skipIfMissingGPU()
    @utx.skipIfMissingFeatures("P3M")
    def test_04_invalid_params_p3m_gpu(self):
        import espressomd.electrostatics

        self.system.time_step = 0.01
        self.add_charged_particles()

        self.check_invalid_params(espressomd.electrostatics.P3MGPU,
                                  mesh=3 * [28], alpha=0.3548, r_cut=4.4434)

    @utx.skipIfMissingFeatures("DP3M")
    def test_04_invalid_params_dp3m_cpu(self):
        import espressomd.magnetostatics

        self.system.time_step = 0.01
        self.add_magnetic_particles()

        self.check_invalid_params(espressomd.magnetostatics.DipolarP3M)

    @utx.skipIfMissingFeatures("P3M")
    def test_04_invalid_params_p3m_elc_cpu(self):
        import espressomd.electrostatics

        self.system.time_step = 0.01
        self.add_charged_particles()

        solver_p3m = espressomd.electrostatics.P3M(
            prefactor=2, accuracy=0.01, tune=False, cao=1,
            r_cut=0.373, alpha=3.81, mesh=(8, 8, 8), check_neutrality=False)
        solver_elc = espressomd.electrostatics.ELC(
            p3m_actor=solver_p3m, gap_size=100, maxPWerror=0.01)
        with self.assertRaisesRegex(ValueError, "gap size too large"):
            self.system.actors.add(solver_elc)

        self.system.actors.clear()
        solver_elc = espressomd.electrostatics.ELC(
            p3m_actor=solver_p3m, gap_size=-1, maxPWerror=0.01)
        with self.assertRaisesRegex(ValueError, "gap_size must be > 0"):
            self.system.actors.add(solver_elc)

        self.system.actors.clear()
        solver_elc = espressomd.electrostatics.ELC(
            p3m_actor=solver_p3m, gap_size=1, maxPWerror=0)
        with self.assertRaisesRegex(ValueError, "maxPWerror must be > 0"):
            self.system.actors.add(solver_elc)

        self.system.part.by_id(0).q = -2

        self.system.actors.clear()
        solver_elc = espressomd.electrostatics.ELC(
            p3m_actor=solver_p3m, gap_size=1,
            maxPWerror=0.01, const_pot=True, pot_diff=3, check_neutrality=False)
        with self.assertRaisesRegex(RuntimeError, "ELC does not currently support non-neutral systems with a dielectric contrast"):
            self.system.actors.add(solver_elc)

        self.system.actors.clear()
        solver_elc = espressomd.electrostatics.ELC(
            p3m_actor=solver_p3m, gap_size=1, maxPWerror=0.01, const_pot=False,
            check_neutrality=False, delta_mid_top=0.5, delta_mid_bot=0.5)
        with self.assertRaisesRegex(RuntimeError, "ELC does not work for non-neutral systems and non-metallic dielectric contrast"):
            self.system.actors.add(solver_elc)

        self.system.part.by_id(0).q = -1

        self.system.actors.clear()
        solver_elc = espressomd.electrostatics.ELC(
            p3m_actor=solver_p3m, gap_size=1, maxPWerror=0.01, const_pot=False,
            check_neutrality=False, delta_mid_top=-1, delta_mid_bot=-1)
        with self.assertRaisesRegex(RuntimeError, "ELC with two parallel metallic boundaries requires the const_pot option"):
            self.system.actors.add(solver_elc)

        self.system.actors.clear()
        solver_dh = espressomd.electrostatics.DH(
            prefactor=1.2, kappa=0.8, r_cut=2.0)
        solver_elc = espressomd.electrostatics.ELC(
            p3m_actor=solver_dh, gap_size=1, maxPWerror=0.01)
        with self.assertRaisesRegex(ValueError, "p3m_actor has to be a P3M solver"):
            self.system.actors.add(solver_elc)

    @utx.skipIfMissingGPU()
    @utx.skipIfMissingFeatures("P3M")
    def test_04_invalid_params_p3m_elc_gpu(self):
        import espressomd.electrostatics

        self.system.time_step = 0.01
        self.add_charged_particles()

        solver_p3m = espressomd.electrostatics.P3MGPU(
            prefactor=2, accuracy=0.01, tune=False, cao=1,
            r_cut=4.4434, alpha=0.3548, mesh=(28, 28, 28))
        solver_elc = espressomd.electrostatics.ELC(
            p3m_actor=solver_p3m, gap_size=1, maxPWerror=0.01)
        with self.assertRaisesRegex(ValueError, "ELC is not set up to work with the GPU P3M"):
            self.system.actors.add(solver_elc)

    @utx.skipIfMissingFeatures("DP3M")
    def test_04_invalid_params_dp3m_dlc_cpu(self):
        import espressomd.magnetostatics
        import espressomd.magnetostatic_extensions

        self.system.time_step = 0.01
        self.add_magnetic_particles()

        dp3m_params = {'accuracy': 1e-6, 'mesh': [25, 25, 25], 'cao': 7,
                       'prefactor': 1.1, 'r_cut': 4.50, 'alpha': 0.8216263}

        solver_dp3m = espressomd.magnetostatics.DipolarP3M(
            epsilon='metallic', tune=False, **dp3m_params)
        self.system.actors.add(solver_dp3m)

        solver_mdlc = espressomd.magnetostatic_extensions.DLC(
            gap_size=100, maxPWerror=1e-5)
        with self.assertRaisesRegex(ValueError, "gap size too large"):
            self.system.actors.add(solver_mdlc)
        self.system.actors.remove(solver_mdlc)

        solver_mdlc = espressomd.magnetostatic_extensions.DLC(
            gap_size=-1, maxPWerror=1e-5)
        with self.assertRaisesRegex(ValueError, "gap_size must be > 0"):
            self.system.actors.add(solver_mdlc)
        self.system.actors.remove(solver_mdlc)

        solver_mdlc = espressomd.magnetostatic_extensions.DLC(
            gap_size=1, maxPWerror=0)
        with self.assertRaisesRegex(ValueError, "maxPWerror must be > 0"):
            self.system.actors.add(solver_mdlc)
        self.system.actors.remove(solver_mdlc)

    ###########################################################
    # block of tests where tuning should not throw exceptions #
    ###########################################################

    def add_actor_assert_failure(self, actor):
        try:
            self.system.actors.add(actor)
        except Exception as err:
            self.fail(f'tuning raised Exception("{err}")')

    @utx.skipIfMissingGPU()
    @utx.skipIfMissingFeatures("P3M")
    def test_09_no_errors_p3m_gpu(self):
        import espressomd.electrostatics

        self.system.time_step = 0.01
        self.add_charged_particles()

        solver = espressomd.electrostatics.P3MGPU(prefactor=2, accuracy=1e-2,
                                                  epsilon='metallic')
        self.add_actor_assert_failure(solver)

    @utx.skipIfMissingFeatures("P3M")
    def test_09_no_errors_p3m_cpu(self):
        import espressomd.electrostatics

        self.system.time_step = 0.01
        self.add_charged_particles()

        solver = espressomd.electrostatics.P3M(prefactor=2, accuracy=0.1)
        valid_params = {
            'mesh_off': solver.default_params()['mesh_off'],  # sentinel
            'cao': 2, 'r_cut': 3.18, 'mesh': 8}

        # tuning with cao or r_cut or mesh constrained, or without constraints
        for key, value in valid_params.items():
            solver = espressomd.electrostatics.P3M(
                prefactor=2, accuracy=1e-2, epsilon=0.0, **{key: value})
            self.add_actor_assert_failure(solver)
            self.system.actors.clear()

    @utx.skipIfMissingFeatures("DP3M")
    def test_09_no_errors_dp3m_cpu(self):
        import espressomd.magnetostatics

        self.system.time_step = 0.01
        self.add_magnetic_particles()

        solver = espressomd.magnetostatics.DipolarP3M(
            prefactor=2, accuracy=0.1)
        valid_params = {
            'mesh_off': solver.default_params()['mesh_off'],  # sentinel
            'cao': 1, 'r_cut': 3.28125, 'mesh': 5}

        # tuning with cao or r_cut or mesh constrained, or without constraints
        for key, value in valid_params.items():
            solver = espressomd.magnetostatics.DipolarP3M(
                prefactor=2, accuracy=1e-2, **{key: value})
            self.add_actor_assert_failure(solver)
            self.system.actors.clear()

    @utx.skipIfMissingFeatures("P3M")
    def test_09_no_errors_p3m_cpu_rescale_mesh(self):
        import espressomd.electrostatics

        self.system.box_l = [10., 15., 20.]
        self.system.time_step = 0.01
        self.add_charged_particles()

        solver = espressomd.electrostatics.P3M(prefactor=2, accuracy=1e-2,
                                               epsilon='metallic',
                                               mesh=[8, -1, -1])
        self.add_actor_assert_failure(solver)
        tuned_mesh = solver.get_params()['mesh']
        self.assertEqual(tuned_mesh[0], 8)
        self.assertEqual(tuned_mesh[1], 12)
        self.assertEqual(tuned_mesh[2], 16)

        # check MD cell reset event
        self.system.box_l = self.system.box_l
        self.system.periodicity = self.system.periodicity

    @utx.skipIfMissingGPU()
    @utx.skipIfMissingFeatures("P3M")
    def test_09_no_errors_p3m_gpu_rescale_mesh(self):
        import espressomd.electrostatics

        self.system.box_l = [10., 15., 20.]
        self.system.time_step = 0.01
        self.add_charged_particles()

        solver = espressomd.electrostatics.P3MGPU(prefactor=2, accuracy=1e-1,
                                                  epsilon='metallic',
                                                  mesh=[20, -1, -1])
        self.add_actor_assert_failure(solver)
        tuned_mesh = solver.get_params()['mesh']
        self.assertEqual(tuned_mesh[0], 20)
        self.assertEqual(tuned_mesh[1], 30)
        self.assertEqual(tuned_mesh[2], 40)

        # check MD cell reset event
        self.system.box_l = self.system.box_l
        self.system.periodicity = self.system.periodicity

    @utx.skipIfMissingFeatures("DP3M")
    def test_09_no_errors_dp3m_cpu_rescale_mesh(self):
        import espressomd.magnetostatics

        self.system.time_step = 0.01
        self.add_magnetic_particles()

        dp3m_params = {'accuracy': 1e-6, 'mesh': [25, 25, 25], 'cao': 7,
                       'prefactor': 1.1, 'r_cut': 4.50, 'alpha': 0.8216263}
        solver_dp3m = espressomd.magnetostatics.DipolarP3M(
            epsilon='metallic', tune=False, **dp3m_params)
        self.add_actor_assert_failure(solver_dp3m)

        # check MD cell reset event
        self.system.box_l = self.system.box_l
        self.system.periodicity = self.system.periodicity


if __name__ == "__main__":
    ut.main()
