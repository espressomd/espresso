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
        self.system.actors.clear()
        self.system.part.clear()
        self.system.box_l = [10., 10., 10.]

    def add_charged_particles(self):
        self.system.part.add(pos=[[0, 0, 0], [.5, .5, .5]], q=[-1, 1])

    def add_magnetic_particles(self):
        self.system.part.add(pos=[[0, 0, 0], [.5, .5, .5]],
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
        with self.assertRaisesRegex(Exception, 'python_p3m_adaptive_tune: ERROR: time_step not set'):
            self.system.actors.add(solver)

    @utx.skipIfMissingFeatures("P3M")
    def test_01_time_not_set_p3m_cpu(self):
        import espressomd.electrostatics

        self.add_charged_particles()

        solver = espressomd.electrostatics.P3M(prefactor=2, accuracy=1e-2)
        with self.assertRaisesRegex(Exception, 'python_p3m_adaptive_tune: ERROR: time_step not set'):
            self.system.actors.add(solver)

    @utx.skipIfMissingFeatures("DP3M")
    def test_01_time_not_set_dp3m_cpu(self):
        import espressomd.magnetostatics

        self.add_magnetic_particles()

        solver = espressomd.magnetostatics.DipolarP3M(
            prefactor=2, accuracy=1e-2)
        with self.assertRaisesRegex(Exception, 'dipolar P3M tuning failed: ERROR: time_step not set'):
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
        with self.assertRaisesRegex(Exception, 'python_p3m_adaptive_tune: ERROR: no charged particles in the system'):
            self.system.actors.add(solver)

    @utx.skipIfMissingFeatures("P3M")
    def test_02_no_particles_p3m_cpu(self):
        import espressomd.electrostatics

        self.system.time_step = 0.01

        solver = espressomd.electrostatics.P3M(prefactor=2, accuracy=1e-2)
        with self.assertRaisesRegex(Exception, 'python_p3m_adaptive_tune: ERROR: no charged particles in the system'):
            self.system.actors.add(solver)

    @utx.skipIfMissingFeatures("DP3M")
    def test_02_no_particles_dp3m_cpu(self):
        import espressomd.magnetostatics

        self.system.time_step = 0.01

        solver = espressomd.magnetostatics.DipolarP3M(
            prefactor=2, accuracy=1e-2)
        with self.assertRaisesRegex(Exception, 'dipolar P3M tuning failed: ERROR: no dipolar particles in the system'):
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
        with self.assertRaisesRegex(Exception, 'python_p3m_adaptive_tune: ERROR: non-metallic epsilon requires cubic box'):
            self.system.actors.add(solver)

    @utx.skipIfMissingFeatures("P3M")
    def test_03_non_cubic_box_p3m_cpu(self):
        import espressomd.electrostatics

        self.system.box_l = [10., 10., 20.]
        self.system.time_step = 0.01
        self.add_charged_particles()

        solver = espressomd.electrostatics.P3M(
            prefactor=2, accuracy=1e-2, epsilon=1)
        with self.assertRaisesRegex(Exception, 'python_p3m_adaptive_tune: ERROR: non-metallic epsilon requires cubic box'):
            self.system.actors.add(solver)

    @utx.skipIfMissingFeatures("DP3M")
    def test_03_non_cubic_box_dp3m_cpu(self):
        import espressomd.magnetostatics

        self.system.box_l = [10., 10., 20.]
        self.system.time_step = 0.01
        self.add_magnetic_particles()

        solver = espressomd.magnetostatics.DipolarP3M(
            prefactor=2, accuracy=1e-2)
        with self.assertRaisesRegex(Exception, 'dipolar P3M tuning failed: ERROR: dipolar P3M requires a cubic box'):
            self.system.actors.add(solver)

    ###########################################################
    # block of tests where tuning should not throw exceptions #
    ###########################################################

    @utx.skipIfMissingGPU()
    @utx.skipIfMissingFeatures("P3M")
    def test_09_no_errors_p3m_gpu(self):
        import espressomd.electrostatics

        self.system.time_step = 0.01
        self.add_charged_particles()

        solver = espressomd.electrostatics.P3MGPU(prefactor=2, accuracy=1e-2,
                                                  epsilon='metallic')
        try:
            self.system.actors.add(solver)
        except Exception as err:
            self.fail('tuning raised Exception("' + str(err) + '")')

    @utx.skipIfMissingFeatures("P3M")
    def test_09_no_errors_p3m_cpu(self):
        import espressomd.electrostatics

        self.system.time_step = 0.01
        self.add_charged_particles()

        solver = espressomd.electrostatics.P3M(prefactor=2, accuracy=1e-2,
                                               epsilon='metallic')
        try:
            self.system.actors.add(solver)
        except Exception as err:
            self.fail('tuning raised Exception("' + str(err) + '")')

    @utx.skipIfMissingFeatures("DP3M")
    def test_09_no_errors_dp3m_cpu(self):
        import espressomd.magnetostatics

        self.system.time_step = 0.01
        self.add_magnetic_particles()

        solver = espressomd.magnetostatics.DipolarP3M(
            prefactor=2, accuracy=1e-2)
        try:
            self.system.actors.add(solver)
        except Exception as err:
            self.fail('tuning raised Exception("' + str(err) + '")')


if __name__ == "__main__":
    ut.main()
