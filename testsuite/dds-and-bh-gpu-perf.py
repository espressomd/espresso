from __future__ import print_function
import time as tm
import unittest as ut
import numpy as np
import espressomd
from espressomd.interactions import *
from espressomd.magnetostatics import *
from espressomd.analyze import *
import math
from tests_common import *
from numpy import linalg as la
from numpy.random import random, seed
from espressomd import assert_features, has_features, missing_features

@ut.skipIf(not has_features(["DIPOLAR_BARNES_HUT"]),
           "Features not available, skipping test!")
class BHGPUPerfTest(ut.TestCase):
    longMessage = True
    # Handle for espresso system
    es = espressomd.System()
    
    def vectorsTheSame(self,a,b):
        tol = 15E-2
        vec_len = la.norm(a - b)
        rel = 2 * vec_len / (la.norm(a) + la.norm(b))
        if rel <= tol:
            return True
        else:
            return False
    
    def stopAll(self):
        for i in range(len(self.es.part)):
             self.es.part[i].v = np.array([0.0,0.0,0.0])
             self.es.part[i].omega_body = np.array([0.0,0.0,0.0])
    
    def run_test_case(self):
        seed(1)
        print("----------------------------------------------")
        print("- Scalability and performance testing: manual durable run only")
        print("- Testcase dds-and-bh-gpu-perf.py")
        print("----------------------------------------------")
        
        pf_bh_gpu = 2.34
        pf_dds_gpu = 3.524
        ratio_dds_gpu_bh_gpu = pf_dds_gpu / pf_bh_gpu
        l = 15
        self.es.box_l = [l, l, l]
        self.es.periodicity = [0, 0, 0]
        self.es.time_step = 1E-4
        self.es.cell_system.skin = 0.1
        
        part_dip = np.zeros((3))
        
        for n in [ 26487, 147543 ]:
            print("{0} particles".format(n))
            force_mag_average = 0.0
            torque_mag_average = 0.0
            dipole_modulus = 1.3
            # scale the box for a large number of particles:
            if n > 1000:
                l *= (n / 541) ** (1 / 3.0)
            for i in range(n):
                part_pos = np.array(random(3)) * l
                costheta = 2 * random() - 1
                sintheta = np.sin(np.arcsin(costheta))
                phi = 2 * np.pi * random()
                part_dip[0] = sintheta * np.cos(phi) * dipole_modulus
                part_dip[1] = sintheta * np.sin(phi) * dipole_modulus
                part_dip[2] = costheta * dipole_modulus
                self.es.part.add(id = i, type = 0, pos = part_pos, dip = part_dip, v = np.array([0,0,0]), omega_body = np.array([0,0,0]))
                
            self.es.non_bonded_inter[0, 0].lennard_jones.set_params(
                epsilon=10.0, sigma=0.5,
                cutoff=0.55, shift="auto")
            self.es.thermostat.set_langevin(kT=0.0, gamma=10.0)
            
            self.es.integrator.set_steepest_descent(f_max=0.0, gamma=0.1, max_displacement=0.1)
            self.es.integrator.run(500)
            self.stopAll()
            self.es.integrator.set_vv()
            
            self.es.non_bonded_inter[0, 0].lennard_jones.set_params(
                epsilon=0.0, sigma=0.0,
                cutoff=-1, shift=0.0)

            self.es.cell_system.skin = 0.0
            self.es.time_step = 0.01
            self.es.thermostat.turn_off()
            
            # gamma should be zero in order to avoid the noise term in force and torque
            self.es.thermostat.set_langevin(kT=1.297, gamma=0.0)
            
            dds_gpu = DipolarDirectSumGpu(prefactor = pf_dds_gpu)
            self.es.actors.add(dds_gpu)
            t1 = tm.time()
            self.es.integrator.run(steps = 0,recalc_forces = True)
            t2 = tm.time()
            dt_dds_gpu = t2 - t1
            
            dds_gpu_f = []
            dds_gpu_t = []
            
            for i in range(n):
                dds_gpu_f.append(self.es.part[i].f)
                dds_gpu_t.append(self.es.part[i].torque_lab)
            dds_gpu_e = Analysis(self.es).energy()["total"]
            
            del dds_gpu
            for i in range(len(self.es.actors.active_actors)):
                self.es.actors.remove(self.es.actors.active_actors[i])
            
            self.es.integrator.run(steps = 0,recalc_forces = True)
            bh_gpu = DipolarBarnesHutGpu(prefactor = pf_bh_gpu, epssq = 400.0, itolsq = 36.0)
            self.es.actors.add(bh_gpu)
            t1 = tm.time()
            self.es.integrator.run(steps = 0,recalc_forces = True)
            t2 = tm.time()
            dt_bh_gpu = t2 - t1
            
            bhgpu_f = []
            bhgpu_t = []
            
            for i in range(n):
                bhgpu_f.append(self.es.part[i].f)
                bhgpu_t.append(self.es.part[i].torque_lab)
            bhgpu_e = Analysis(self.es).energy()["total"]
            
            for i in range(n):
                force_mag_average += la.norm(dds_gpu_f[i])
                torque_mag_average += la.norm(dds_gpu_t[i])
            
            force_mag_average /= n
            torque_mag_average /= n
            
            cutoff = 1E-2
            
            # compare
            for i in range(n):
                if la.norm(dds_gpu_t[i]) > cutoff * torque_mag_average:
                    self.assertTrue(self.vectorsTheSame(np.array(dds_gpu_t[i]),ratio_dds_gpu_bh_gpu * np.array(bhgpu_t[i])), \
                                    msg = 'Torques on particle do not match. i={0} dds_gpu_t={1} ratio_dds_gpu_bh_gpu*bhgpu_t={2}'.format(i,np.array(dds_gpu_t[i]), ratio_dds_gpu_bh_gpu * np.array(bhgpu_t[i])))
                if la.norm(dds_gpu_f[i]) > cutoff * force_mag_average:
                    self.assertTrue(self.vectorsTheSame(np.array(dds_gpu_f[i]),ratio_dds_gpu_bh_gpu * np.array(bhgpu_f[i])), \
                                    msg = 'Forces on particle do not match: i={0} dds_gpu_f={1} ratio_dds_gpu_bh_gpu*bhgpu_f={2}'.format(i,np.array(dds_gpu_f[i]), ratio_dds_gpu_bh_gpu * np.array(bhgpu_f[i])))
            self.assertTrue(abs(dds_gpu_e - bhgpu_e * ratio_dds_gpu_bh_gpu) <= abs(1E-3 * dds_gpu_e), \
                            msg = 'Energies for dawaanr {0} and dds_gpu {1} do not match.'.format(dds_gpu_e,ratio_dds_gpu_bh_gpu * bhgpu_e))
            
            print("=== Performance comparison ===")
            print("dt_dds_gpu = {0}".format(dt_dds_gpu))
            print("dt_bh_gpu = {0}".format(dt_bh_gpu))
            
            self.es.integrator.run(steps = 0,recalc_forces = True)
            
            del bh_gpu
            for i in range(len(self.es.actors.active_actors)):
                self.es.actors.remove(self.es.actors.active_actors[i])
            self.es.part.clear()

    def test(self):
        if (self.es.cell_system.get_state()["n_nodes"] > 1):
            print("NOTE: Ignoring testcase for n_nodes > 1")
        else:
            self.run_test_case()
    
if __name__ == '__main__':
    print("Features: ", espressomd.features())
    ut.main()
