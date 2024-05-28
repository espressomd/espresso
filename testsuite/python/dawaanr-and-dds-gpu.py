#
# Copyright (C) 2010-2024 The ESPResSo project
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
import tests_common
import numpy as np

import espressomd
import espressomd.magnetostatics


@utx.skipIfMissingGPU()
@utx.skipIfMissingFeatures(["DIPOLAR_DIRECT_SUM", "LENNARD_JONES"])
class DDSGPUTest(ut.TestCase):
    system = espressomd.System(box_l=[1.0, 1.0, 1.0])
    np.random.seed(seed=42)

    def test(self):
        pf_dds_cpu = 2.34
        pf_dds_gpu = 3.524
        ratio = pf_dds_cpu / pf_dds_gpu
        system = self.system
        system.box_l = [15., 15., 15.]
        system.periodicity = [False, False, False]
        system.time_step = 0.01
        system.cell_system.skin = 0.1

        for n in [128, 541]:
            dipole_modulus = 1.3
            part_dip = dipole_modulus * tests_common.random_dipoles(n)
            part_pos = np.random.random((n, 3)) * system.box_l[0]
            system.part.add(pos=part_pos, dip=part_dip)

            # warmup
            system.non_bonded_inter[0, 0].lennard_jones.set_params(
                epsilon=10.0, sigma=0.5, cutoff=0.55, shift="auto")
            system.thermostat.set_langevin(kT=0.0, gamma=10.0, seed=42)
            system.integrator.set_vv()
            system.integrator.run(steps=20)
            system.non_bonded_inter[0, 0].lennard_jones.deactivate()
            system.thermostat.turn_off()
            system.galilei.kill_particle_motion(rotation=True)

            # gamma should be zero in order to avoid the noise term in force
            # and torque
            system.thermostat.set_langevin(kT=1.297, gamma=0.0)

            dds_cpu = espressomd.magnetostatics.DipolarDirectSumCpu(
                prefactor=pf_dds_cpu)
            system.magnetostatics.solver = dds_cpu
            system.integrator.run(steps=0, recalc_forces=True)

            dawaanr_f = np.copy(system.part.all().f)
            dawaanr_t = np.copy(system.part.all().torque_lab)
            dawaanr_e = system.analysis.energy()["total"]

            del dds_cpu
            system.magnetostatics.clear()

            system.integrator.run(steps=0, recalc_forces=True)
            dds_gpu = espressomd.magnetostatics.DipolarDirectSumGpu(
                prefactor=pf_dds_gpu)
            system.magnetostatics.solver = dds_gpu
            system.integrator.run(steps=0, recalc_forces=True)

            dds_gpu_f = np.copy(system.part.all().f)
            dds_gpu_t = np.copy(system.part.all().torque_lab)
            dds_gpu_e = system.analysis.energy()["total"]

            np.testing.assert_allclose(dawaanr_t, ratio * dds_gpu_t,
                                       atol=1e-10, rtol=1e-2)
            np.testing.assert_allclose(dawaanr_f, ratio * dds_gpu_f,
                                       atol=1e-10, rtol=1e-2)
            np.testing.assert_allclose(dawaanr_e, ratio * dds_gpu_e,
                                       atol=1e-10, rtol=1e-5)

            # check MD cell reset has no impact
            system.change_volume_and_rescale_particles(system.box_l[0], "x")
            system.periodicity = system.periodicity
            system.cell_system.node_grid = system.cell_system.node_grid
            system.integrator.run(steps=0, recalc_forces=True)

            del dds_gpu
            system.magnetostatics.clear()
            system.part.clear()


if __name__ == '__main__':
    ut.main()
