from __future__ import print_function
import unittest as ut
import numpy as np
from numpy.random import random, seed
import espressomd
import math


@ut.skipIf(not espressomd.has_features(["MASS",
                                        "PARTICLE_ANISOTROPY",
                                        "ROTATIONAL_INERTIA",
                                        "LANGEVIN_PER_PARTICLE"]),
           "Features not available, skipping test!")
class ThermoTest(ut.TestCase):
    longMessage = True
    # Handle for espresso system
    es = espressomd.System()

    def run_test_case(self, test_case):
        seed(1)
        # Decelleration / drift tests
        self.es.time_step = 0.007
        # gamma_tran/gamma_rot matrix: [2 types of particless] x [3 dimensions
        # X Y Z]
        gamma_tran = np.zeros((2, 3))
        gamma_tr = np.zeros((2, 3))
        gamma_rot = np.zeros((2, 3))
        gamma_rot_validate = np.zeros((2, 3))
        # Global gamma for tests without particle-specific gammas:
        gamma_global = np.ones((3))
        # Per-paricle values for the remaining tests:
        # Either translational friction isotropy is required
        # or both translational and rotational ones.
        # Otherwise these types of motion will interfere.
        # ..Let's test both cases depending on the particle index.
        gamma_tran[0, 0] = np.array(0.5 + random(1))
        gamma_tran[0, 1] = gamma_tran[0, 0]
        gamma_tran[0, 2] = gamma_tran[0, 0]
        if test_case < 4:
            gamma_rot[0, :] = np.array((0.5 + random(3)) * 2.0 / 3.0)
        else:
            if "BROWNIAN_DYNAMICS" in espressomd.features():
                # Isotropy is required here for drift tests (see below)
                gamma_rot[0, 0] = np.array((0.5 + random(1)) * 2.0 / 3.0)
                gamma_rot[0, 1] = gamma_rot[0, 0]
                gamma_rot[0, 2] = gamma_rot[0, 0]
        
        gamma_tran[1, 0] = np.array(0.5 + random(1))
        gamma_tran[1, 1] = gamma_tran[1, 0]
        gamma_tran[1, 2] = gamma_tran[1, 0]
        gamma_rot[1, 0] = np.array((0.5 + random(1)) * 2.0 / 3.0)
        gamma_rot[1, 1] = gamma_rot[1, 0]
        gamma_rot[1, 2] = gamma_rot[1, 0]

        if test_case < 4:
            self.es.thermostat.set_langevin(
                kT=0.0,
                gamma=[
                    gamma_global[0],
                    gamma_global[1],
                    gamma_global[2]])
        else:
            if "BROWNIAN_DYNAMICS" in espressomd.features():
                # Brownian thermostat activation
                self.es.thermostat.turn_off()
                self.es.thermostat.set_brownian(
                    kT=0.0,
                    gamma=[
                        gamma_global[0],
                        gamma_global[1],
                        gamma_global[2]])

        self.es.cell_system.skin = 5.0
        mass = 12.74
        J = [10.0, 10.0, 10.0]

        for i in range(len(self.es.part)):
            self.es.part[i].remove()

        for i in range(2):
            self.es.part.add(rotation=(1,1,1), pos=np.array([0.0, 0.0, 0.0]), id=i)
            self.es.part[i].v = np.array([1.0, 1.0, 1.0])
            if "ROTATION" in espressomd.features():
                self.es.part[i].omega_body = np.array([1.0, 1.0, 1.0])
            self.es.part[i].mass = mass
            self.es.part[i].rinertia = np.array(J)

        print("\n")

        if test_case == 0:
            print("================================================")
            print("Group of test cases 0-3: Langevin thermostat")
            print("================================================")
            print("------------------------------------------------")
            print("Test " + str(test_case) + ": no particle specific values")
            print("------------------------------------------------")
            # No assignments are needed.

        if test_case == 1:
            print("------------------------------------------------")
            print("Test " + str(test_case) +
                  ": particle specific gamma but not temperature")
            print("------------------------------------------------")
            if "PARTICLE_ANISOTROPY" in espressomd.features():
                self.es.part[0].gamma = gamma_tran[0, :]
                self.es.part[1].gamma = gamma_tran[1, :]
            else:
                self.es.part[0].gamma = gamma_tran[0, 0]
                self.es.part[1].gamma = gamma_tran[1, 0]
            if "ROTATION" in espressomd.features():
                self.es.part[0].gamma_rot = gamma_rot[0, :]
                self.es.part[1].gamma_rot = gamma_rot[1, :]

        if test_case == 2:
            print("------------------------------------------------")
            print("Test " + str(test_case) +
                  ": particle specific temperature but not gamma")
            print("------------------------------------------------")
            self.es.part[0].temp = 0.0
            self.es.part[1].temp = 0.0

        if test_case == 3:
            print("------------------------------------------------")
            print("Test " + str(test_case) +
                  ": both particle specific gamma and temperature")
            print("------------------------------------------------")
            self.es.part[0].temp = 0.0
            self.es.part[1].temp = 0.0
            if "PARTICLE_ANISOTROPY" in espressomd.features():
                self.es.part[0].gamma = gamma_tran[0, :]
                self.es.part[1].gamma = gamma_tran[1, :]
            else:
                self.es.part[0].gamma = gamma_tran[0, 0]
                self.es.part[1].gamma = gamma_tran[1, 0]
            if "ROTATION" in espressomd.features():
                self.es.part[0].gamma_rot = gamma_rot[0, :]
                self.es.part[1].gamma_rot = gamma_rot[1, :]

        if "BROWNIAN_DYNAMICS" in espressomd.features():
            if test_case == 4:
                print("================================================")
                print("Group of test cases 4-7: Brownian thermostat")
                print("================================================")
                print("------------------------------------------------")
                print("Test " + str(test_case) + ": no particle specific values")
                print("------------------------------------------------")
                # No assignments are needed.
    
            if test_case == 5:
                print("------------------------------------------------")
                print("Test " + str(test_case) +
                    ": particle specific gamma but not temperature")
                print("------------------------------------------------")
                if "PARTICLE_ANISOTROPY" in espressomd.features():
                    self.es.part[0].gamma = gamma_tran[0, :]
                    self.es.part[1].gamma = gamma_tran[1, :]
                else:
                    self.es.part[0].gamma = gamma_tran[0, 0]
                    self.es.part[1].gamma = gamma_tran[1, 0]
                if "ROTATION" in espressomd.features():
                    self.es.part[0].gamma_rot = gamma_rot[0, :]
                    self.es.part[1].gamma_rot = gamma_rot[1, :]
    
            if test_case == 6:
                print("------------------------------------------------")
                print("Test " + str(test_case) +
                    ": particle specific temperature but not gamma")
                print("------------------------------------------------")
                self.es.part[0].temp = 0.0
                self.es.part[1].temp = 0.0
    
            if test_case == 7:
                print("------------------------------------------------")
                print("Test " + str(test_case) +
                    ": both particle specific gamma and temperature")
                print("------------------------------------------------")
                self.es.part[0].temp = 0.0
                self.es.part[1].temp = 0.0
                if "PARTICLE_ANISOTROPY" in espressomd.features():
                    self.es.part[0].gamma = gamma_tran[0, :]
                    self.es.part[1].gamma = gamma_tran[1, :]
                else:
                    self.es.part[0].gamma = gamma_tran[0, 0]
                    self.es.part[1].gamma = gamma_tran[1, 0]
                if "ROTATION" in espressomd.features():
                    self.es.part[0].gamma_rot = gamma_rot[0, :]
                    self.es.part[1].gamma_rot = gamma_rot[1, :]

        if test_case == 1 or test_case == 3 or test_case == 5 or test_case == 7:
            gamma_tr = gamma_tran
            gamma_rot_validate = gamma_rot
        else:
            for k in range(2):
                gamma_tr[k, :] = gamma_global[:]
                gamma_rot_validate[k, :] = gamma_global[:]

        self.es.time = 0.0

        if test_case < 4:
            self.es.time_step = 0.007
            tol = 1.25E-4
            for i in range(100):
                for k in range(3):
                    self.assertLess(
                        abs(self.es.part[0].v[k] - math.exp(- gamma_tr[0, k] * self.es.time / mass)), tol)
                    self.assertLess(
                        abs(self.es.part[1].v[k] - math.exp(- gamma_tr[1, k] * self.es.time / mass)), tol)
                    if "ROTATION" in espressomd.features():
                        self.assertLess(abs(
                            self.es.part[0].omega_body[k] - math.exp(- gamma_rot_validate[0, k] * self.es.time / J[k])), tol)
                        self.assertLess(abs(
                            self.es.part[1].omega_body[k] - math.exp(- gamma_rot_validate[1, k] * self.es.time / J[k])), tol)
                self.es.integrator.run(10)
        # Brownian dynamics drift test
        else:
            tol = 1.25E-4
            # Brownian dynamics is designed for larger time-steps
            # propagate for time_step >> max(tau0_tran, tau0_rot)
            # tau0_tran = mass / gamma = 12.74 / 1.0
            # tau0_rot = I / gamma = 10 / 1.0
            # Hence, let's select time_step ~ 100
            self.es.time_step = 100
            if "EXTERNAL_FORCES" in espressomd.features():
                f0 = np.array([-1.2, 58.3578, 0.002])
                f1 = np.array([-15.112, -2.0, 368])
                self.es.part[0].ext_force = f0
                self.es.part[1].ext_force = f1
                if "ROTATION" in espressomd.features():
                    tor0 = np.array([12, 0.022, 87])
                    tor1 = np.array([-0.03, -174, 368])
                    self.es.part[0].ext_torque = tor0
                    self.es.part[1].ext_torque = tor1
                    # Let's set the dipole perpendicular to the torque
                    if "DIPOLES" in espressomd.features():
                        dip0 = np.array([0.0, tor0[2], -tor0[1]])
                        dip1 = np.array([-tor1[2], 0.0, tor1[0]])
                        self.es.part[0].dip = dip0
                        self.es.part[1].dip = dip1
                for i in range(100):
                    self.es.integrator.run(10)
                    #print("i = {0}".format(i))
                    for k in range(3):
                        self.assertLess(
                            abs(self.es.part[0].v[k] - f0[k] / gamma_tr[0, k]), tol)
                        self.assertLess(
                            abs(self.es.part[1].v[k] - f1[k] / gamma_tr[1, k]), tol)
                        self.assertLess(
                            abs(self.es.part[0].pos[k] - self.es.time * f0[k] / gamma_tr[0, k]), tol)
                        self.assertLess(
                            abs(self.es.part[1].pos[k] - self.es.time * f1[k] / gamma_tr[1, k]), tol)
                        if "ROTATION" in espressomd.features():
                            self.assertLess(abs(
                                self.es.part[0].omega_lab[k] - tor0[k] / gamma_rot_validate[0, k]), tol)
                            self.assertLess(abs(
                                self.es.part[1].omega_lab[k] - tor1[k] / gamma_rot_validate[1, k]), tol)
                    if "ROTATION" in espressomd.features() and "DIPOLES" in espressomd.features():
                        cos_alpha0 = np.dot(dip0,self.es.part[0].dip) / (np.linalg.norm(dip0) * self.es.part[0].dipm)
                        cos_alpha0_test = np.cos(self.es.time * np.linalg.norm(tor0) / gamma_rot_validate[0, 0])
                        
                        cos_alpha1 = np.dot(dip1,self.es.part[1].dip) / (np.linalg.norm(dip1) * self.es.part[1].dipm)
                        cos_alpha1_test = np.cos(self.es.time * np.linalg.norm(tor1) / gamma_rot_validate[1, 0])
                        
                        #print("cos_alpha0 = {0}, cos_alpha0_test={1}".format(cos_alpha0, cos_alpha0_test))
                        #print("dip0 = {0}, self.es.part[0].dip={1}".format(dip0, self.es.part[0].dip))
                        self.assertLess(abs(cos_alpha0 - cos_alpha0_test), tol)
                        self.assertLess(abs(cos_alpha1 - cos_alpha1_test), tol)

        for i in range(len(self.es.part)):
            self.es.part[i].remove()

        # thermalization
        # Checks if every degree of freedom has 1/2 kT of energy, even when
        # mass and inertia tensor are active

        # 2 different langevin parameters for particles
        temp = np.array([2.5, 2.0])
        D_tr = np.zeros((2, 3))
        for k in range(2):
            if "PARTICLE_ANISOTROPY" in espressomd.features():
                gamma_tran[k, :] = np.array((0.4 + random(3)) * 10)
                gamma_rot[k, :] = np.array((0.2 + random(3)) * 20)
            else:
                gamma_tran[k, 0] = np.array((0.4 + random(1)) * 10)
                gamma_tran[k, 1] = gamma_tran[k, 0]
                gamma_tran[k, 2] = gamma_tran[k, 0]
                gamma_rot[k, 0] = np.array((0.2 + random(1)) * 20)
                gamma_rot[k, 1] = gamma_rot[k, 0]
                gamma_rot[k, 2] = gamma_rot[k, 0]

        box = 10.0
        self.es.box_l = [box, box, box]
        if espressomd.has_features(("PARTIAL_PERIODIC",)):
            self.es.periodicity = 0, 0, 0
        kT = 1.5
        gamma_global = np.ones((3))

        if test_case in [2,3,6,7]:
            halfkT = temp / 2.0
        else:
            halfkT = np.array([kT, kT]) / 2.0

        if test_case in [1,3,5,7]:
            gamma_tr = gamma_tran
        else:
            for k in range(2):
                gamma_tr[k, :] = gamma_global[:]

        # translational diffusion
        for k in range(2):
            D_tr[k, :] = 2.0 * halfkT[k] / gamma_tr[k, :]

        if "PARTICLE_ANISOTROPY" in espressomd.features():
            if test_case < 4:
                self.es.thermostat.set_langevin(
                    kT=kT,
                    gamma=[
                        gamma_global[0],
                        gamma_global[1],
                        gamma_global[2]])
            else:
                if "BROWNIAN_DYNAMICS" in espressomd.features():
                    self.es.thermostat.turn_off()
                    self.es.thermostat.set_brownian(
                        kT=kT,
                        gamma=[
                            gamma_global[0],
                            gamma_global[1],
                            gamma_global[2]])
        else:
            if test_case < 4:
                self.es.thermostat.set_langevin(kT=kT, gamma=gamma_global[0])
            else:
                if "BROWNIAN_DYNAMICS" in espressomd.features():
                    self.es.thermostat.turn_off()
                    self.es.thermostat.set_brownian(kT=kT, gamma=gamma_global[0])

        # no need to rebuild Verlet lists, avoid it
        self.es.cell_system.skin = 5.0
        if test_case < 4:
            self.es.time_step = 0.03
        else:
            self.es.time_step = 10
        n = 200
        mass = (0.2 + random()) * 7.0
        J = np.array((0.2 + random(3)) * 7.0)

        for i in range(n):
            for k in range(2):
                ind = i + k * n
                part_pos = np.array(random(3) * box)
                part_v = np.array([0.0, 0.0, 0.0])
                part_omega_body = np.array([0.0, 0.0, 0.0])
                self.es.part.add(rotation=(1,1,1), id=ind, mass=mass, rinertia=J,
                                 pos=part_pos, v=part_v)
                if "ROTATION" in espressomd.features():
                    self.es.part[ind].omega_body = part_omega_body
                if test_case in [1,5]:
                    if "PARTICLE_ANISOTROPY" in espressomd.features():
                        self.es.part[ind].gamma = gamma_tran[k, :]
                    else:
                        self.es.part[ind].gamma = gamma_tran[k, 0]
                    if "ROTATION" in espressomd.features():
                        self.es.part[ind].gamma_rot = gamma_rot[k, :]

                if test_case in [2,6]:
                    self.es.part[ind].temp = temp[k]
                if test_case in [3,7]:
                    if "PARTICLE_ANISOTROPY" in espressomd.features():
                        self.es.part[ind].gamma = gamma_tran[k, :]
                    else:
                        self.es.part[ind].gamma = gamma_tran[k, 0]
                    if "ROTATION" in espressomd.features():
                        self.es.part[ind].gamma_rot = gamma_rot[k, :]
                    self.es.part[ind].temp = temp[k]

        # Get rid of short range calculations if exclusions are on
        # if espressomd.has_features("EXCLUSIONS"):

        # matrices: [2 types of particless] x [3 dimensions X Y Z]
        # velocity^2, omega^2, position^2
        v2 = np.zeros((2, 3))
        o2 = np.zeros((2, 3))
        dr2 = np.zeros((2, 3))
        sigma2_tr = np.zeros((2))
        dr_norm = np.zeros((2))

        pos0 = np.zeros((2 * n, 3))
        for p in range(n):
            for k in range(2):
                ind = p + k * n
                pos0[ind, :] = self.es.part[ind].pos
        dt0 = mass / gamma_tr

        loops = 200
        print("Thermalizing...")
        therm_steps = 150
        self.es.integrator.run(therm_steps)
        print("Measuring...")

        int_steps = 5
        for i in range(loops):
            self.es.integrator.run(int_steps)
            # Get kinetic energy in each degree of freedom for all particles
            for p in range(n):
                for k in range(2):
                    ind = p + k * n
                    v = self.es.part[ind].v
                    if "ROTATION" in espressomd.features():
                        o = self.es.part[ind].omega_body
                        o2[k, :] = o2[k, :] + np.power(o[:], 2)
                    pos = self.es.part[ind].pos
                    v2[k, :] = v2[k, :] + np.power(v[:], 2)
                    dr2[k, :] = np.power((pos[:] - pos0[ind, :]), 2)
                    dt = (int_steps * (i + 1) + therm_steps) * \
                        self.es.time_step
                    # translational diffusion variance: after a closed-form
                    # integration of the Langevin EOM
                    sigma2_tr[k] = 0.0
                    for j in range(3):
                        sigma2_tr[k] = sigma2_tr[k] + D_tr[k,
                                                           j] * (2 * dt + dt0[k,
                                                                              j] * (- 3 + 4 * math.exp(- dt / dt0[k,
                                                                                                                  j]) - math.exp(- 2 * dt / dt0[k,
                                                                                                                                                j])))
                    dr_norm[k] = dr_norm[k] + \
                        (sum(dr2[k, :]) - sigma2_tr[k]) / sigma2_tr[k]

        tolerance = 0.15
        Ev = 0.5 * mass * v2 / (n * loops)
        Eo = 0.5 * J * o2 / (n * loops)
        dv = np.zeros((2))
        do = np.zeros((2))
        do_vec = np.zeros((2, 3))
        for k in range(2):
            dv[k] = sum(Ev[k, :]) / (3 * halfkT[k]) - 1.0
            do[k] = sum(Eo[k, :]) / (3 * halfkT[k]) - 1.0
            do_vec[k, :] = Eo[k, :] / halfkT[k] - 1.0
        dr_norm = dr_norm / (n * loops)
        for k in range(2):
            print("\n")
            print("k = " + str(k))
            print("mass = " + str(mass))
            print("gamma_tr = {0} {1} {2}".format(
                gamma_tr[k, 0], gamma_tr[k, 1], gamma_tr[k, 2]))
            if test_case in [1,5,3,7]:
                print("gamma_rot = {0} {1} {2}".format(
                    gamma_rot[k, 0], gamma_rot[k, 1], gamma_rot[k, 2]))
            else:
                print("gamma_global = {0} {1} {2}".format(
                    gamma_global[0], gamma_global[1], gamma_global[2]))
            print("Moment of inertia principal components: = " + str(J))
            print("1/2 kT = " + str(halfkT[k]))
            print("Translational energy: {0} {1} {2}".format(
                Ev[k, 0], Ev[k, 1], Ev[k, 2]))
            print("Rotational energy: {0} {1} {2}".format(
                Eo[k, 0], Eo[k, 1], Eo[k, 2]))

            print("Deviation in translational energy: " + str(dv[k]))
            if "ROTATION" in espressomd.features():
                print("Deviation in rotational energy: " + str(do[k]))
                print("Deviation in rotational energy per degrees of freedom: {0} {1} {2}".format(
                    do_vec[k, 0], do_vec[k, 1], do_vec[k, 2]))
            print(
                "Deviation in translational diffusion: {0} ".format(
                    dr_norm[k]))

            self.assertLessEqual(
                abs(
                    dv[k]),
                tolerance,
                msg='Relative deviation in translational energy too large: {0}'.format(
                    dv[k]))
            if "ROTATION" in espressomd.features():
                self.assertLessEqual(
                    abs(
                        do[k]),
                    tolerance,
                    msg='Relative deviation in rotational energy too large: {0}'.format(
                        do[k]))
                self.assertLessEqual(abs(
                    do_vec[k, 0]), tolerance, msg='Relative deviation in rotational energy per the body axis X is too large: {0}'.format(do_vec[k, 0]))
                self.assertLessEqual(abs(
                    do_vec[k, 1]), tolerance, msg='Relative deviation in rotational energy per the body axis Y is too large: {0}'.format(do_vec[k, 1]))
                self.assertLessEqual(abs(
                    do_vec[k, 2]), tolerance, msg='Relative deviation in rotational energy per the body axis Z is too large: {0}'.format(do_vec[k, 2]))
            self.assertLessEqual(
                abs(
                    dr_norm[k]),
                tolerance,
                msg='Relative deviation in translational diffusion is too large: {0}'.format(
                    dr_norm[k]))

    def test(self):
        if not "BROWNIAN_DYNAMICS" in espressomd.features():
            n_test_cases = 4
        else:
            n_test_cases = 8
        for i in range(n_test_cases):
            self.run_test_case(i)


if __name__ == '__main__':
    print("Features: ", espressomd.features())
    ut.main()
