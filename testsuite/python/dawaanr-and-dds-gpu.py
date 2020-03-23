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
from numpy.random import random
import numpy as np

import espressomd
import espressomd.interactions
import espressomd.magnetostatics
import espressomd.analyze


@utx.skipIfMissingGPU()
@utx.skipIfMissingFeatures(["DIPOLES", "ROTATION", "LENNARD_JONES"])
class DDSGPUTest(ut.TestCase):
    # Handle for espresso system
    es = espressomd.System(box_l=[1.0, 1.0, 1.0])

    def stopAll(self):
        for i in range(len(self.es.part)):
            self.es.part[i].v = np.array([0.0, 0.0, 0.0])
            self.es.part[i].omega_body = np.array([0.0, 0.0, 0.0])

    @ut.skipIf(es.cell_system.get_state()["n_nodes"] > 1,
               "Skipping test: only runs for n_nodes == 1")
    def test(self):
        pf_dds_gpu = 2.34
        pf_dawaanr = 3.524
        ratio_dawaanr_dds_gpu = pf_dawaanr / pf_dds_gpu
        l = 15
        self.es.box_l = [l, l, l]
        self.es.periodicity = [0, 0, 0]
        self.es.time_step = 1E-4
        self.es.cell_system.skin = 0.1

        part_dip = np.zeros((3))

        for n in [128, 541]:
            dipole_modulus = 1.3
            for i in range(n):
                part_pos = np.array(random(3)) * l
                costheta = 2 * random() - 1
                sintheta = np.sin(np.arcsin(costheta))
                phi = 2 * np.pi * random()
                part_dip[0] = sintheta * np.cos(phi) * dipole_modulus
                part_dip[1] = sintheta * np.sin(phi) * dipole_modulus
                part_dip[2] = costheta * dipole_modulus
                self.es.part.add(id=i, type=0, pos=part_pos, dip=part_dip,
                                 v=np.array([0, 0, 0]), omega_body=np.array([0, 0, 0]))

            self.es.non_bonded_inter[0, 0].lennard_jones.set_params(
                epsilon=10.0, sigma=0.5, cutoff=0.55, shift="auto")
            self.es.thermostat.set_langevin(kT=0.0, gamma=10.0, seed=42)

            self.es.integrator.set_steepest_descent(
                f_max=0.0, gamma=0.1, max_displacement=0.1)
            self.es.integrator.run(500)
            self.stopAll()
            self.es.integrator.set_vv()

            self.es.non_bonded_inter[0, 0].lennard_jones.set_params(
                epsilon=0.0, sigma=0.0, cutoff=0.0, shift=0.0)

            self.es.cell_system.skin = 0.0
            self.es.time_step = 0.01
            self.es.thermostat.turn_off()
            # gamma should be zero in order to avoid the noise term in force
            # and torque
            self.es.thermostat.set_langevin(kT=1.297, gamma=0.0, seed=42)

            dds_cpu = espressomd.magnetostatics.DipolarDirectSumCpu(
                prefactor=pf_dawaanr)
            self.es.actors.add(dds_cpu)
            self.es.integrator.run(steps=0, recalc_forces=True)

            dawaanr_f = []
            dawaanr_t = []

            for i in range(n):
                dawaanr_f.append(self.es.part[i].f)
                dawaanr_t.append(self.es.part[i].torque_lab)
            dawaanr_e = self.es.analysis.energy()["total"]

            del dds_cpu
            for i in range(len(self.es.actors.active_actors)):
                self.es.actors.remove(self.es.actors.active_actors[i])

            self.es.integrator.run(steps=0, recalc_forces=True)
            dds_gpu = espressomd.magnetostatics.DipolarDirectSumGpu(
                prefactor=pf_dds_gpu)
            self.es.actors.add(dds_gpu)
            self.es.integrator.run(steps=0, recalc_forces=True)

            ddsgpu_f = []
            ddsgpu_t = []

            for i in range(n):
                ddsgpu_f.append(self.es.part[i].f)
                ddsgpu_t.append(self.es.part[i].torque_lab)
            ddsgpu_e = self.es.analysis.energy()["total"]

            # compare
            for i in range(n):
                np.testing.assert_allclose(
                    np.array(dawaanr_t[i]),
                    ratio_dawaanr_dds_gpu * np.array(ddsgpu_t[i]),
                    err_msg='Torques on particle do not match for particle {}'
                    .format(i), atol=3e-3)
                np.testing.assert_allclose(
                    np.array(dawaanr_f[i]),
                    ratio_dawaanr_dds_gpu * np.array(ddsgpu_f[i]),
                    err_msg='Forces on particle do not match for particle i={}'
                    .format(i), atol=3e-3)
            self.assertAlmostEqual(
                dawaanr_e,
                ddsgpu_e * ratio_dawaanr_dds_gpu,
                places=2,
                msg='Energies for dawaanr {0} and dds_gpu {1} do not match.'
                .format(dawaanr_e, ratio_dawaanr_dds_gpu * ddsgpu_e))

            self.es.integrator.run(steps=0, recalc_forces=True)

            del dds_gpu
            self.es.actors.clear()
            self.es.part.clear()


if __name__ == '__main__':
    ut.main()
