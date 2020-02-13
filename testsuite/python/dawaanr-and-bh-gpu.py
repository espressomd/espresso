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
import numpy as np

import espressomd
import espressomd.magnetostatics
import espressomd.analyze
import espressomd.cuda_init


def stopAll(system):
    system.part[:].v = np.zeros(3)
    system.part[:].omega_body = np.zeros(3)


@utx.skipIfMissingGPU()
@utx.skipIfMissingFeatures(["DIPOLAR_BARNES_HUT", "LENNARD_JONES"])
class BHGPUTest(ut.TestCase):
    system = espressomd.System(box_l=[1, 1, 1])
    np.random.seed(seed=42)

    def vectorsTheSame(self, a, b):
        tol = 5E-2
        vec_len = np.linalg.norm(a - b)
        rel = 2 * vec_len / (np.linalg.norm(a) + np.linalg.norm(b))
        return rel <= tol

    @ut.skipIf(system.cell_system.get_state()["n_nodes"] > 1,
               "Skipping test: only runs for n_nodes == 1")
    def test(self):
        np.random.seed(1)
        pf_bh_gpu = 2.34
        pf_dawaanr = 3.524
        ratio_dawaanr_bh_gpu = pf_dawaanr / pf_bh_gpu
        l = 15
        self.system.box_l = [l, l, l]
        self.system.periodicity = [0, 0, 0]
        self.system.time_step = 1E-4
        self.system.cell_system.skin = 0.1

        part_dip = np.zeros((3))

        for n in [128, 541]:
            dipole_modulus = 1.3
            # scale the box for a large number of particles:
            if n > 1000:
                l *= (n / 541) ** (1 / 3.0)
            for i in range(n):
                part_pos = np.array(np.random.random(3)) * l
                costheta = 2 * np.random.random() - 1
                sintheta = np.sin(np.arcsin(costheta))
                phi = 2 * np.pi * np.random.random()
                part_dip[0] = sintheta * np.cos(phi) * dipole_modulus
                part_dip[1] = sintheta * np.sin(phi) * dipole_modulus
                part_dip[2] = costheta * dipole_modulus
                self.system.part.add(id=i, type=0, pos=part_pos, dip=part_dip,
                                     v=np.array([0, 0, 0]), omega_body=np.array([0, 0, 0]))

            self.system.non_bonded_inter[0, 0].lennard_jones.set_params(
                epsilon=10.0, sigma=0.5, cutoff=0.55, shift="auto")
            self.system.thermostat.set_langevin(kT=0.0, gamma=10.0, seed=42)
            stopAll(self.system)
            self.system.integrator.set_vv()

            self.system.non_bonded_inter[0, 0].lennard_jones.set_params(
                epsilon=0.0, sigma=0.0, cutoff=-1, shift=0.0)

            self.system.cell_system.skin = 0.0
            self.system.time_step = 0.01
            self.system.thermostat.turn_off()

            # gamma should be zero in order to avoid the noise term in force
            # and torque
            self.system.thermostat.set_langevin(kT=1.297, gamma=0.0)

            dds_cpu = espressomd.magnetostatics.DipolarDirectSumCpu(
                prefactor=pf_dawaanr)
            self.system.actors.add(dds_cpu)
            self.system.integrator.run(steps=0, recalc_forces=True)

            dawaanr_f = []
            dawaanr_t = []

            for i in range(n):
                dawaanr_f.append(self.system.part[i].f)
                dawaanr_t.append(self.system.part[i].torque_lab)
            dawaanr_e = espressomd.analyze.Analysis(
                self.system).energy()["total"]

            del dds_cpu
            self.system.actors.clear()

            self.system.integrator.run(steps=0, recalc_forces=True)
            bh_gpu = espressomd.magnetostatics.DipolarBarnesHutGpu(
                prefactor=pf_bh_gpu, epssq=200.0, itolsq=8.0)
            self.system.actors.add(bh_gpu)
            self.system.integrator.run(steps=0, recalc_forces=True)

            bhgpu_f = []
            bhgpu_t = []

            for i in range(n):
                bhgpu_f.append(self.system.part[i].f)
                bhgpu_t.append(self.system.part[i].torque_lab)
            bhgpu_e = espressomd.analyze.Analysis(
                self.system).energy()["total"]

            # compare
            for i in range(n):
                self.assertTrue(
                    self.vectorsTheSame(
                        np.array(dawaanr_t[i]),
                        ratio_dawaanr_bh_gpu * np.array(bhgpu_t[i])),
                    msg='Torques on particle do not match. i={0} dawaanr_t={1} '
                        'ratio_dawaanr_bh_gpu*bhgpu_t={2}'.format(
                            i, np.array(dawaanr_t[i]),
                            ratio_dawaanr_bh_gpu * np.array(bhgpu_t[i])))
                self.assertTrue(
                    self.vectorsTheSame(
                        np.array(dawaanr_f[i]),
                        ratio_dawaanr_bh_gpu * np.array(bhgpu_f[i])),
                    msg='Forces on particle do not match: i={0} dawaanr_f={1} '
                        'ratio_dawaanr_bh_gpu*bhgpu_f={2}'.format(
                            i, np.array(dawaanr_f[i]),
                            ratio_dawaanr_bh_gpu * np.array(bhgpu_f[i])))
            self.assertLessEqual(
                abs(dawaanr_e - bhgpu_e * ratio_dawaanr_bh_gpu),
                abs(1E-3 * dawaanr_e),
                msg='Energies for dawaanr {0} and bh_gpu {1} do not match.'
                .format(dawaanr_e, ratio_dawaanr_bh_gpu * bhgpu_e))

            self.system.integrator.run(steps=0, recalc_forces=True)

            del bh_gpu
            self.system.actors.clear()
            self.system.part.clear()


if __name__ == '__main__':
    ut.main()
