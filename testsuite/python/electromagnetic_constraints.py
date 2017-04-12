from __future__ import print_function
import unittest as ut
import numpy as np
from numpy.random import random
import espressomd
import math
from espressomd.interactions import *
from espressomd.magnetostatics import *
from espressomd.electrostatics import *
from espressomd.shapes import *
from espressomd.constraints import Constraint

if "MASS" in espressomd.features() and "ROTATIONAL_INERTIA" in espressomd.features() \
and "DIPOLES" in espressomd.features() and "PARTIAL_PERIODIC" in espressomd.features() \
and "CONSTRAINTS" in espressomd.features() and "ELECTROSTATICS" in espressomd.features() :
    class ESMSDriftTest(ut.TestCase):
        longMessage = True
        # Handle for espresso system
        es = espressomd.System()
        n_nodes = 0

        def run_test_case(self):
            box = 1000
            self.es.box_l = [box, box, box]
            self.es.time_step = 7E-3
            self.es.cell_system.set_layered(n_layers = int(math.ceil(25.0 / self.n_nodes)))
            n = 9
            gamma0 = 7.984
            kT = 2.64
            self.es.periodicity = [1, 1, 0]
            self.es.cell_system.skin = 0
            self.es.thermostat.set_langevin(kT = kT, gamma = gamma0)
            J = np.array([4.384, 4.384, 4.384])
            mass = 3.847
            B_mag = 10.0
            q = 1.0
            sigma = 5.0 / (2.0 * np.pi)
            
            # This structure has the center of mass zc
            zc = box / 2.0 + 7 / 3.0
            self.es.part.add(pos = np.array([box / 2.0, box / 2.0, box / 2.0]),id = 0)
            self.es.part.add(pos = np.array([box / 2.0, box / 2.0, box / 2.0 + 2.0]),id = 1)
            self.es.part.add(pos = np.array([box / 2.0, box / 2.0, box / 2.0 + 4.0]),id = 2)
            self.es.part.add(pos = np.array([box / 2.0, box / 2.0 + 2.0, box / 2.0 + 1.0]),id = 3)
            self.es.part.add(pos = np.array([box / 2.0, box / 2.0 + 2.0, box / 2.0 + 3.0]),id = 4)
            self.es.part.add(pos = np.array([box / 2.0, box / 2.0 + 2.0, box / 2.0 + 5.0]),id = 5)
            self.es.part.add(pos = np.array([box / 2.0, box / 2.0 + 4.0, box / 2.0]),id = 6)
            self.es.part.add(pos = np.array([box / 2.0, box / 2.0 + 4.0, box / 2.0 + 2.0]),id = 7)
            self.es.part.add(pos = np.array([box / 2.0, box / 2.0 + 4.0, box / 2.0 + 4.0]),id = 8)
            
            for p in range(n):
                theta = np.pi * random()
                phi = 2 * np.pi * random()
                self.es.part[p].type = 0
                self.es.part[p].mass = mass
                self.es.part[p].q = q
                self.es.part[p].dip = np.array([0.0, 0.0, 1.0])
                self.es.part[p].rinertia = J
                self.es.part[p].v = np.zeros((3))
                self.es.part[p].omega_body = np.zeros((3))

            #Interactions
            dds_cpu = DipolarDirectSumCpu(bjerrum_length = 20.0 / kT)
            self.es.actors.add(dds_cpu)
            mmm2d = MMM2D(bjerrum_length = 1.0 / kT, maxPWerror = 1E-6)
            self.es.actors.add(mmm2d)
            wall1 = Wall(normal = [0.0, 0.0, 1.0], dist = 0.02 * box)
            constraint1 = Constraint(shape = wall1, only_positive = 0, particle_type = 2, penetrable = 1, ext_electric_field = 0.0, ext_magn_field = B_mag)
            wall2 = Wall(normal = [0.0, 0.0, 1.0], dist = 0.03 * box)
            constraint2 = Constraint(shape = wall2, only_positive = 0, particle_type = 2, penetrable = 0, ext_electric_field = 0.0, ext_magn_field = 0.0)
            wall3 = Wall(normal = [0.0, 0.0, 1.0], dist = 0.01 * box)
            constraint3 = Constraint(shape = wall3, only_positive = 0, particle_type = 2, penetrable = 1, ext_electric_field = 2.0 * np.pi * sigma, ext_magn_field = 0.0)
            self.es.constraints.add(constraint1)
            self.es.constraints.add(constraint2)
            self.es.constraints.add(constraint3)
            
            sig = 5.0 / 6.0
            cut = 1.4 * 1.12246 * sig
            eps = 5
            shift = 0.25 * eps
            self.es.non_bonded_inter[0, 0].lennard_jones.set_params(
                epsilon=eps, sigma=sig,
                cutoff=cut, shift=0.25 * eps)

            sig = 5.0 / 6.0
            cut = 0.4 * 1.12246 * sig
            eps = 5
            shift = 0.25 * eps
            self.es.non_bonded_inter[0, 2].lennard_jones.set_params(
                epsilon=eps, sigma=sig,
                cutoff=cut, shift=0.25 * eps)

            q_tot = q * n
            F = 2 * np.pi * sigma * q_tot
            gamma_tot = gamma0 * n
            mass_tot = mass * n
            
            self.es.time = 0.0
            
            for i in range(600):
                self.es.integrator.run(10)

            # The aggregate macroscopic parameters
            # Center of mass Z component
            pos_c_z = 0

            for p in range(n):
                pos_c_z += self.es.part[p].pos[2]
            
            t = self.es.time
            pos_c_z = pos_c_z / n - box / 2.0
            pos_c_z_exp = (F /  gamma_tot**2) * (mass_tot * (math.exp(- gamma_tot * t / mass_tot) - 1) + gamma_tot * t) + 7 / 3.0
            
            tol = 0.2
            self.assertTrue(abs(pos_c_z - pos_c_z_exp) / pos_c_z_exp <= tol, msg = 'Z-drift deviation is too large: {0} vs the expected value {1}'.format(pos_c_z,pos_c_z_exp))

        def test(self):
            self.n_nodes = self.es.cell_system.get_state()["n_nodes"]
            print("----------------------------------------------------------------")
            print("- Testcase esms-driven-aggregate-drift.py running on {0} nodes -".format(self.n_nodes))
            print("----------------------------------------------------------------")
            self.run_test_case()

if __name__ == '__main__':
    print("Features: ", espressomd.features())
    ut.main()
