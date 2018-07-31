from __future__ import print_function
import unittest as ut
import numpy as np
import espressomd
from espressomd.interactions import *
from espressomd.magnetostatics import *
from espressomd.analyze import *
import math
from tests_common import *
from numpy import linalg as la
from numpy.random import random
from espressomd import has_features


@ut.skipIf(not has_features(["DIPOLES",
                             "CUDA",
                             "PARTIAL_PERIODIC",
                             "ROTATION"]),
           "Features not available, skipping test!")
class DDSGPUTest(ut.TestCase):
    longMessage = True
    # Handle for espresso system
    es = espressomd.System(box_l=[1.0, 1.0, 1.0])
    es.seed  = es.cell_system.get_state()['n_nodes'] * [1234]


    def stopAll(self):
        for i in range(len(self.es.part)):
            self.es.part[i].v = np.array([0.0, 0.0, 0.0])
            self.es.part[i].omega_body = np.array([0.0, 0.0, 0.0])

    def run_test_case(self):
        print("----------------------------------------------")
        print("- Testcase dawaanr-and-dds-gpu.py")
        print("----------------------------------------------")

        pf_dds_gpu = 2.34
        pf_dawaanr = 3.524
        ratio_dawaanr_dds_gpu = pf_dawaanr / pf_dds_gpu
        l = 15
        self.es.box_l = [l, l, l]
        self.es.periodicity = [0, 0, 0]
        self.es.time_step = 1E-4
        self.es.cell_system.skin = 0.1

        part_dip = np.zeros((3))

        for n in [110, 111, 540, 541]:
            print("{0} particles".format(n))
            dipole_modulus = 1.3
            for i in range(n):
                part_pos = np.array(random(3)) * l
                costheta = 2 * random() - 1
                sintheta = np.sin(np.arcsin(costheta))
                phi = 2 * np.pi * random()
                part_dip[0] = sintheta * np.cos(phi) * dipole_modulus
                part_dip[1] = sintheta * np.sin(phi) * dipole_modulus
                part_dip[2] = costheta * dipole_modulus
                self.es.part.add(id=i, type=0, pos=part_pos, dip=part_dip, v=np.array(
                    [0, 0, 0]), omega_body=np.array([0, 0, 0]))

            self.es.non_bonded_inter[0, 0].lennard_jones.set_params(
                epsilon=10.0, sigma=0.5,
                cutoff=0.55, shift="auto")
            self.es.thermostat.set_langevin(kT=0.0, gamma=10.0)

            self.es.integrator.set_steepest_descent(
                f_max=0.0, gamma=0.1, max_displacement=0.1)
            self.es.integrator.run(500)
            self.stopAll()
            self.es.integrator.set_vv()

            self.es.non_bonded_inter[0, 0].lennard_jones.set_params(
                epsilon=0.0, sigma=0.0,
                cutoff=0.0, shift=0.0)

            self.es.cell_system.skin = 0.0
            self.es.time_step = 0.01
            self.es.thermostat.turn_off()
            # gamma should be zero in order to avoid the noise term in force
            # and torque
            self.es.thermostat.set_langevin(kT=1.297, gamma=0.0)

            dds_cpu = DipolarDirectSumCpu(prefactor=pf_dawaanr)
            self.es.actors.add(dds_cpu)
            self.es.integrator.run(steps=0, recalc_forces=True)

            dawaanr_f = []
            dawaanr_t = []

            for i in range(n):
                dawaanr_f.append(self.es.part[i].f)
                dawaanr_t.append(self.es.part[i].torque_lab)
            dawaanr_e = Analysis(self.es).energy()["total"]

            del dds_cpu
            for i in range(len(self.es.actors.active_actors)):
                self.es.actors.remove(self.es.actors.active_actors[i])

            self.es.integrator.run(steps=0, recalc_forces=True)
            dds_gpu = DipolarDirectSumGpu(prefactor=pf_dds_gpu)
            self.es.actors.add(dds_gpu)
            self.es.integrator.run(steps=0, recalc_forces=True)

            ddsgpu_f = []
            ddsgpu_t = []

            for i in range(n):
                ddsgpu_f.append(self.es.part[i].f)
                ddsgpu_t.append(self.es.part[i].torque_lab)
            ddsgpu_e = Analysis(self.es).energy()["total"]

            # compare
            for i in range(n):
                np.testing.assert_allclose(np.array(dawaanr_t[i]),
                        ratio_dawaanr_dds_gpu * np.array(ddsgpu_t[i]),
                    err_msg='Torques on particle do not match for particle {}'.format(i), atol=3e-3)
                np.testing.assert_allclose(np.array(dawaanr_f[i]),
                        ratio_dawaanr_dds_gpu * np.array(ddsgpu_f[i]),
                    err_msg='Forces on particle do not match for particle i={}'.format(i), atol=3e-3)
            self.assertAlmostEqual(
                dawaanr_e,
                ddsgpu_e *
                ratio_dawaanr_dds_gpu,
                places=2,
                msg='Energies for dawaanr {0} and dds_gpu {1} do not match.'.format(
                    dawaanr_e,
                    ratio_dawaanr_dds_gpu *
                    ddsgpu_e))

            self.es.integrator.run(steps=0, recalc_forces=True)

            del dds_gpu
            for i in range(len(self.es.actors.active_actors)):
                self.es.actors.remove(self.es.actors.active_actors[i])
            # for i in reversed(range(len(self.es.part))):
            #    self.es.part[i].remove()
            self.es.part.clear()

    def test(self):
        if (self.es.cell_system.get_state()["n_nodes"] > 1):
            print("NOTE: Ignoring testcase for n_nodes > 1")
        else:
            self.run_test_case()


if __name__ == '__main__':
    print("Features: ", espressomd.features())
    ut.main()
