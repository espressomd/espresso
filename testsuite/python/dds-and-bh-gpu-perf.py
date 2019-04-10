# Copyright (C) 2010-2018 The ESPResSo project
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
from __future__ import print_function

import time as tm
import unittest as ut
import math
import numpy as np
from numpy import linalg as la
from numpy.random import random, seed

import espressomd
import espressomd.magnetostatics


@ut.skipIf(
    not espressomd.gpu_available() or not espressomd.has_features(
        ["DIPOLAR_BARNES_HUT"]),
           "Features or gpu not available, skipping test!")
class BHGPUPerfTest(ut.TestCase):
    # Handle for espresso system
    system = espressomd.System(box_l=[10.0, 10.0, 10.0])

    def vectorsTheSame(self, a, b):
        tol = 15E-2
        vec_len = la.norm(a - b)
        rel = 2 * vec_len / (la.norm(a) + la.norm(b))
        return rel <= tol

    def stopAll(self):
        for i in range(len(self.system.part)):
            self.system.part[i].v = np.array([0.0, 0.0, 0.0])
            self.system.part[i].omega_body = np.array([0.0, 0.0, 0.0])

    def run_test_case(self):
        seed(1)
    
        pf_bh_gpu = 2.34
        pf_dds_gpu = 3.524
        ratio_dds_gpu_bh_gpu = pf_dds_gpu / pf_bh_gpu
        l = 15
        self.system.periodicity = [0, 0, 0]
        self.system.time_step = 1E-4
        self.system.cell_system.skin = 0.1

        part_dip = np.zeros((3))

        for n in [26487, 147543]:
            force_mag_average = 0.0
            torque_mag_average = 0.0
            dipole_modulus = 1.3
            # scale the box for a large number of particles:
            if n > 1000:
                l *= 1.5 * (n / 541) ** (1 / 3.0)
            self.system.box_l = [l, l, l]
            for i in range(n):
                part_pos = np.array(random(3)) * l
                costheta = 2 * random() - 1
                sintheta = np.sin(np.arcsin(costheta))
                phi = 2 * np.pi * random()
                part_dip[0] = sintheta * np.cos(phi) * dipole_modulus
                part_dip[1] = sintheta * np.sin(phi) * dipole_modulus
                part_dip[2] = costheta * dipole_modulus
                self.system.part.add(id=i, type=0, pos=part_pos, dip=part_dip, v=np.array(
                    [0, 0, 0]), omega_body=np.array([0, 0, 0]))

            self.system.non_bonded_inter[0, 0].lennard_jones.set_params(
                epsilon=10.0, sigma=0.5,
                cutoff=0.55, shift="auto")
            self.system.thermostat.set_langevin(kT=0.0, gamma=10.0, seed=42)

            self.system.integrator.set_steepest_descent(
                f_max=0.0, gamma=0.1, max_displacement=0.1)
            self.system.integrator.run(500)
            self.stopAll()
            self.system.integrator.set_vv()

            self.system.non_bonded_inter[0, 0].lennard_jones.set_params(
                epsilon=0.0, sigma=0.0,
                cutoff=-1, shift=0.0)

            self.system.cell_system.skin = 0.0
            self.system.time_step = 0.01
            self.system.thermostat.turn_off()

            # gamma should be zero in order to avoid the noise term in force
            # and torque
            self.system.thermostat.set_langevin(kT=1.297, gamma=0.0, seed=42)

            dds_gpu = espressomd.magnetostatics.DipolarDirectSumGpu(
                prefactor=pf_dds_gpu)
            self.system.actors.add(dds_gpu)
            t1 = tm.time()
            self.system.integrator.run(steps=0, recalc_forces=True)
            t2 = tm.time()
            dt_dds_gpu = t2 - t1

            dds_gpu_f = []
            dds_gpu_t = []

            for i in range(n):
                dds_gpu_f.append(self.system.part[i].f)
                dds_gpu_t.append(self.system.part[i].torque_lab)
            dds_gpu_e = self.system.analysis.energy()["total"]

            del dds_gpu
            self.system.actors.clear()

            self.system.integrator.run(steps=0, recalc_forces=True)
            bh_gpu = espressomd.magnetostatics.DipolarBarnesHutGpu(
                prefactor=pf_bh_gpu, epssq=400.0, itolsq=36.0)
            self.system.actors.add(bh_gpu)
            t1 = tm.time()
            self.system.integrator.run(steps=0, recalc_forces=True)
            t2 = tm.time()
            dt_bh_gpu = t2 - t1

            bhgpu_f = []
            bhgpu_t = []

            for i in range(n):
                bhgpu_f.append(self.system.part[i].f)
                bhgpu_t.append(self.system.part[i].torque_lab)
            bhgpu_e = self.system.analysis.energy()["total"]

            for i in range(n):
                force_mag_average += la.norm(dds_gpu_f[i])
                torque_mag_average += la.norm(dds_gpu_t[i])

            force_mag_average /= n
            torque_mag_average /= n

            cutoff = 1E-2

            # compare
            for i in range(n):
                if la.norm(dds_gpu_t[i]) > cutoff * torque_mag_average:
                    self.assertTrue(
                        self.vectorsTheSame(
                            np.array(
                                dds_gpu_t[i]), ratio_dds_gpu_bh_gpu * np.array(
                                    bhgpu_t[i])),
                                    msg='Torques on particle do not match. i={0} dds_gpu_t={1} ratio_dds_gpu_bh_gpu*bhgpu_t={2}'.format(i, np.array(dds_gpu_t[i]), ratio_dds_gpu_bh_gpu * np.array(bhgpu_t[i])))
                if la.norm(dds_gpu_f[i]) > cutoff * force_mag_average:
                    self.assertTrue(
                        self.vectorsTheSame(
                            np.array(
                                dds_gpu_f[i]), ratio_dds_gpu_bh_gpu * np.array(
                                    bhgpu_f[i])),
                                    msg='Forces on particle do not match: i={0} dds_gpu_f={1} ratio_dds_gpu_bh_gpu*bhgpu_f={2}'.format(i, np.array(dds_gpu_f[i]), ratio_dds_gpu_bh_gpu * np.array(bhgpu_f[i])))
            self.assertTrue(
                abs(dds_gpu_e - bhgpu_e * ratio_dds_gpu_bh_gpu) <= abs(
                    1E-3 * dds_gpu_e),
                            msg='Energies for dawaanr {0} and dds_gpu {1} do not match.'.format(dds_gpu_e, ratio_dds_gpu_bh_gpu * bhgpu_e))

            print("=== Performance comparison ===")
            print("dt_dds_gpu = {0}".format(dt_dds_gpu))
            print("dt_bh_gpu = {0}".format(dt_bh_gpu))

            self.system.integrator.run(steps=0, recalc_forces=True)

            del bh_gpu
            for i in range(len(self.system.actors.active_actors)):
                self.system.actors.remove(self.system.actors.active_actors[i])
            self.system.part.clear()

    def test(self):
        if (self.system.cell_system.get_state()["n_nodes"] > 1):
            print("NOTE: Ignoring testcase for n_nodes > 1")
        else:
            self.run_test_case()

if __name__ == '__main__':
    ut.main()
