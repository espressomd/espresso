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
import unittest as ut
import unittest_decorators as utx
import tests_common
import numpy as np

import espressomd
import espressomd.interactions
import espressomd.magnetostatics
import espressomd.analyze
import espressomd.galilei


@utx.skipIfMissingGPU()
@utx.skipIfMissingFeatures(["DIPOLES", "ROTATION", "LENNARD_JONES"])
class DDSGPUTest(ut.TestCase):
    # Handle for espresso system
    system = espressomd.System(box_l=[1.0, 1.0, 1.0])

    @ut.skipIf(system.cell_system.get_state()["n_nodes"] > 1,
               "Skipping test: only runs for n_nodes == 1")
    def test(self):
        pf_dds_gpu = 2.34
        pf_dawaanr = 3.524
        ratio_dawaanr_dds_gpu = pf_dawaanr / pf_dds_gpu
        self.system.box_l = 3 * [15]
        self.system.periodicity = [0, 0, 0]
        self.system.time_step = 1E-4
        self.system.cell_system.skin = 0.1

        for n in [128, 541]:
            dipole_modulus = 1.3
            part_dip = dipole_modulus * tests_common.random_dipoles(n)
            part_pos = np.random.random((n, 3)) * self.system.box_l[0]
            self.system.part.add(pos=part_pos, dip=part_dip)

            self.system.non_bonded_inter[0, 0].lennard_jones.set_params(
                epsilon=10.0, sigma=0.5, cutoff=0.55, shift="auto")

            self.system.thermostat.turn_off()
            self.system.integrator.set_steepest_descent(
                f_max=0.0, gamma=0.1, max_displacement=0.1)
            self.system.integrator.run(500)
            g = espressomd.galilei.GalileiTransform()
            g.kill_particle_motion(rotation=True)
            self.system.integrator.set_vv()

            self.system.non_bonded_inter[0, 0].lennard_jones.set_params(
                epsilon=0.0, sigma=0.0, cutoff=0.0, shift=0.0)

            self.system.cell_system.skin = 0.0
            self.system.time_step = 0.01
            self.system.thermostat.turn_off()
            # gamma should be zero in order to avoid the noise term in force
            # and torque
            self.system.thermostat.set_langevin(kT=1.297, gamma=0.0, seed=42)

            dds_cpu = espressomd.magnetostatics.DipolarDirectSumCpu(
                prefactor=pf_dawaanr)
            self.system.actors.add(dds_cpu)
            self.system.integrator.run(steps=0, recalc_forces=True)

            dawaanr_f = np.copy(self.system.part.all().f)
            dawaanr_t = np.copy(self.system.part.all().torque_lab)
            dawaanr_e = self.system.analysis.energy()["total"]

            del dds_cpu
            for i in range(len(self.system.actors.active_actors)):
                self.system.actors.remove(self.system.actors.active_actors[i])

            self.system.integrator.run(steps=0, recalc_forces=True)
            dds_gpu = espressomd.magnetostatics.DipolarDirectSumGpu(
                prefactor=pf_dds_gpu)
            self.system.actors.add(dds_gpu)
            self.system.integrator.run(steps=0, recalc_forces=True)

            ddsgpu_f = np.copy(self.system.part.all().f)
            ddsgpu_t = np.copy(self.system.part.all().torque_lab)
            ddsgpu_e = self.system.analysis.energy()["total"]

            # compare
            for i in range(n):
                np.testing.assert_allclose(
                    np.array(dawaanr_t[i]),
                    ratio_dawaanr_dds_gpu * np.array(ddsgpu_t[i]),
                    err_msg=f'Torques do not match for particle {i}',
                    atol=3e-3)
                np.testing.assert_allclose(
                    np.array(dawaanr_f[i]),
                    ratio_dawaanr_dds_gpu * np.array(ddsgpu_f[i]),
                    err_msg=f'Forces do not match for particle {i}',
                    atol=3e-3)
            self.assertAlmostEqual(
                dawaanr_e,
                ddsgpu_e * ratio_dawaanr_dds_gpu,
                places=2,
                msg='Energies for dawaanr {0} and dds_gpu {1} do not match.'
                .format(dawaanr_e, ratio_dawaanr_dds_gpu * ddsgpu_e))

            self.system.integrator.run(steps=0, recalc_forces=True)

            del dds_gpu
            self.system.actors.clear()
            self.system.part.clear()


if __name__ == '__main__':
    ut.main()
