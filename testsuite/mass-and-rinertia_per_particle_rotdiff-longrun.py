from __future__ import print_function
import unittest as ut
import numpy as np
from numpy.random import random, seed
import espressomd
import math
import tests_common as tc

@ut.skipIf(not espressomd.has_features(["MASS",
                                        "PARTICLE_ANISOTROPY",
                                        "ROTATIONAL_INERTIA",
                                        "LANGEVIN_PER_PARTICLE"]),
           "Features not available, skipping test!")
class ThermoTest(ut.TestCase):
    longMessage = True
    # Handle for espresso system
    es = espressomd.System(box_l=[1.0E5, 1.0E5, 1.0E5])
    rot_flag = 0
    
    def set_thermo_all(self, test_case, kT, gamma_global, gamma_global_rot):
        if test_case < 4:
            self.es.thermostat.set_langevin(
                kT=kT,
                gamma=[
                    gamma_global[0],
                    gamma_global[1],
                    gamma_global[2]])
        elif test_case == 4 and self.rot_flag == 1:
            self.es.thermostat.set_langevin(
                kT=kT,
                gamma=[
                    gamma_global[0],
                    gamma_global[1],
                    gamma_global[2]],
                gamma_rotation=[
                    gamma_global_rot[0],
                    gamma_global_rot[1],
                    gamma_global_rot[2]])
        elif test_case < 8 + self.rot_flag:
            if "BROWNIAN_DYNAMICS" in espressomd.features():
                # Brownian thermostat activation
                self.es.thermostat.turn_off()
                self.es.thermostat.set_brownian(
                    kT=kT,
                    gamma=[
                        gamma_global[0],
                        gamma_global[1],
                        gamma_global[2]])
        elif test_case == 9:
            if "BROWNIAN_DYNAMICS" in espressomd.features():
                # Brownian thermostat activation
                self.es.thermostat.turn_off()
                self.es.thermostat.set_brownian(
                    kT=kT,
                    gamma=[
                        gamma_global[0],
                        gamma_global[1],
                        gamma_global[2]],
                    gamma_rotation=[
                        gamma_global_rot[0],
                        gamma_global_rot[1],
                        gamma_global_rot[2]])

    def run_test_case(self, test_case):
        seed(2)
        # Decelleration
        self.es.time_step = 0.007
        self.es.part.clear()
        # gamma_tran/gamma_rot matrix: [2 types of particless] x [3 dimensions
        # X Y Z]
        gamma_tran = np.zeros((2, 3))
        gamma_tr = np.zeros((2, 3))
        gamma_rot = np.zeros((2, 3))
        gamma_rot_validate = np.zeros((2, 3))
        # The translational gamma isotropy is required here. See comments below.
        # Global gamma for tests without particle-specific gammas:
        # gamma_global = np.ones((3))
        gamma_global = np.zeros((3))
        gamma_global[0] = np.array((0.5 + np.random.random()) * 2.0 / 3.0)
        gamma_global[1] = gamma_global[0]
        gamma_global[2] = gamma_global[0]
        # Additional test case (5th) for the specific global rotational gamma set.
        gamma_global_rot = np.array((0.5 + np.random.random(3)) * 2.0 / 3.0)
        # Per-paricle values for the remaining decelleration tests:
        # Either translational friction isotropy is required
        # or both translational and rotational ones.
        # Otherwise these types of motion will interfere.
        # ..Let's test both cases depending on the particle index.
        gamma_tran[0, 0] = np.array(0.5 + np.random.random())
        gamma_tran[0, 1] = gamma_tran[0, 0]
        gamma_tran[0, 2] = gamma_tran[0, 0]
        if test_case < 4 + self.rot_flag:
            gamma_rot[0, :] = np.array((0.5 + random(3)) * 2.0 / 3.0)
        else:
            if "BROWNIAN_DYNAMICS" in espressomd.features():
                # Isotropy is required here for the drag tests (see below)
                gamma_rot[0, 0] = np.array((0.5 + random(1)) * 2.0 / 3.0)
                gamma_rot[0, 1] = gamma_rot[0, 0]
                gamma_rot[0, 2] = gamma_rot[0, 0]
        
        gamma_tran[1, 0] = np.array(0.5 + np.random.random())
        gamma_tran[1, 1] = gamma_tran[1, 0]
        gamma_tran[1, 2] = gamma_tran[1, 0]
        gamma_rot[1, 0] = np.array((0.5 + np.random.random()) * 2.0 / 3.0)
        gamma_rot[1, 1] = gamma_rot[1, 0]
        gamma_rot[1, 2] = gamma_rot[1, 0]
        
        self.set_thermo_all(test_case, 0.0, gamma_global, gamma_global_rot)

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

        for k in range(2):
            if test_case == 0:
                if (k == 0):
                    print("================================================")
                    print("Langevin thermostat test cases")
                    print("================================================")
                    print("------------------------------------------------")
                    print("Test " + str(test_case) + ": no particle specific values")
                    print("------------------------------------------------")
                # No assignments are needed.
    
            if test_case == 1:
                if (k == 0):
                    print("------------------------------------------------")
                    print("Test " + str(test_case) +
                        ": particle specific gamma but not temperature")
                    print("------------------------------------------------")
                self.es.part[k].gamma = gamma_tran[k, :]
                if "ROTATION" in espressomd.features():
                    self.es.part[k].gamma_rot = gamma_rot[k, :]
    
            if test_case == 2:
                if (k == 0):
                    print("------------------------------------------------")
                    print("Test " + str(test_case) +
                        ": particle specific temperature but not gamma")
                    print("------------------------------------------------")
                self.es.part[k].temp = 0.0
    
            if test_case == 3:
                if (k == 0):
                    print("------------------------------------------------")
                    print("Test " + str(test_case) +
                        ": both particle specific gamma and temperature")
                    print("------------------------------------------------")
                self.es.part[k].temp = 0.0
                self.es.part[k].gamma = gamma_tran[k, :]
                if "ROTATION" in espressomd.features():
                    self.es.part[k].gamma_rot = gamma_rot[k, :]
            
            if test_case == 4 and self.rot_flag == 1:
                if (k == 0):
                    print("------------------------------------------------")
                    print("Test " + str(test_case) + ": no particle specific values.")
                    print("Rotational specific global thermostat")
                    print("------------------------------------------------")
                    # No assignments are needed.

        if "BROWNIAN_DYNAMICS" in espressomd.features():
            for k in range(2):
                if test_case == 4 + self.rot_flag:
                    if (k == 0):
                        print("================================================")
                        print("Brownian thermostat test cases")
                        print("================================================")
                        print("------------------------------------------------")
                        print("Test " + str(test_case) + ": no particle specific values")
                        print("------------------------------------------------")
                    # No assignments are needed.
        
                if test_case == 5 + self.rot_flag:
                    if (k == 0):
                        print("------------------------------------------------")
                        print("Test " + str(test_case) +
                            ": particle specific gamma but not temperature")
                        print("------------------------------------------------")
                    if "PARTICLE_ANISOTROPY" in espressomd.features():
                        self.es.part[k].gamma = gamma_tran[k, :]
                    else:
                        self.es.part[k].gamma = gamma_tran[k, 0]
                    if "ROTATION" in espressomd.features():
                        self.es.part[k].gamma_rot = gamma_rot[k, :]
        
                if test_case == 6 + self.rot_flag:
                    if (k == 0):
                        print("------------------------------------------------")
                        print("Test " + str(test_case) +
                            ": particle specific temperature but not gamma")
                        print("------------------------------------------------")
                    self.es.part[k].temp = 0.0
        
                if test_case == 7 + self.rot_flag:
                    if (k == 0):
                        print("------------------------------------------------")
                        print("Test " + str(test_case) +
                            ": both particle specific gamma and temperature")
                        print("------------------------------------------------")
                    self.es.part[k].temp = 0.0
                    if "PARTICLE_ANISOTROPY" in espressomd.features():
                        self.es.part[k].gamma = gamma_tran[k, :]
                    else:
                        self.es.part[k].gamma = gamma_tran[k, 0]
                    if "ROTATION" in espressomd.features():
                        self.es.part[k].gamma_rot = gamma_rot[k, :]
                        
                if test_case == 8 + self.rot_flag:
                    if (k == 0):
                        if (k == 0):
                            print("------------------------------------------------")
                            print("Test " + str(test_case) + ": no particle specific values.")
                            print("Rotational specific global thermostat")
                            print("------------------------------------------------")
                        # No assignments are needed.

        if test_case == 1 or test_case == 3 or test_case == (5 + self.rot_flag) or test_case == (7 + self.rot_flag):
            gamma_tr = gamma_tran
            gamma_rot_validate = gamma_rot
        else:
            for k in range(2):
                gamma_tr[k, :] = gamma_global[:]
                if (test_case == 4 or test_case == 9) and self.rot_flag == 1:
                    gamma_rot_validate[k, :] = gamma_global_rot[:]
                else:
                    gamma_rot_validate[k, :] = gamma_global[:]

        if test_case < 4 + self.rot_flag:
        # Langevin thermostat only. Brownian thermostat is defined
        # over larger time-step by its concept.
            self.es.time = 0.0
            tol = 1.25E-4
            for i in range(100):
                for k in range(2):
                    for j in range(3):
                        # Note: velocity is defined in the lab frame of reference
                        # while gamma_tr is defined in the body one.
                        # Hence, only isotropic gamma_tr could be tested here.
                        self.assertLess(
                            abs(self.es.part[k].v[j] - math.exp(- gamma_tr[k, j] * self.es.time / mass)), tol)
                        if "ROTATION" in espressomd.features():
                            self.assertLess(abs(
                                self.es.part[k].omega_body[j] - math.exp(- gamma_rot_validate[k, j] * self.es.time / J[j])), tol)
                self.es.integrator.run(10)
        # The drag terminal velocity tests
        ##################################
        #..aka the drift in case of the electrostatics
        
        # Isotropic reassignment is required here for the drag tests
        gamma_global_rot = np.zeros((3))
        gamma_global_rot[0] = np.array((0.5 + np.random.random()) * 2.0 / 3.0)
        gamma_global_rot[1] = gamma_global_rot[0]
        gamma_global_rot[2] = gamma_global_rot[0]
        # Second particle already has isotropic gamma_rot
        # A correction is needed for the 1st one:
        gamma_rot[0, 0] = np.array((0.5 + random(1)) * 2.0 / 3.0)
        gamma_rot[0, 1] = gamma_rot[0, 0]
        gamma_rot[0, 2] = gamma_rot[0, 0]
        
        if test_case == 1 or test_case == 3 or test_case == (5 + self.rot_flag) or test_case == (7 + self.rot_flag):
            gamma_tr = gamma_tran
            gamma_rot_validate = gamma_rot
            # A correction is needed for the 1st particle only:
            if "ROTATION" in espressomd.features():
                self.es.part[0].gamma_rot = gamma_rot[0, :]
        else:
            for k in range(2):
                gamma_tr[k, :] = gamma_global[:]
                if (test_case == 4 or test_case == 9) and self.rot_flag == 1:
                    gamma_rot_validate[k, :] = gamma_global_rot[:]
                else:
                    gamma_rot_validate[k, :] = gamma_global[:]

        self.set_thermo_all(test_case, 0.0, gamma_global, gamma_global_rot)

        self.es.time = 0.0
        self.es.time_step = 7E-5
        # The terminal velocity is starting since t >> t0 = mass / gamma
        t0_max = -1.0
        for k in range(2):
            t0_max = max(t0_max, max(mass / gamma_tr[k, :]), max(J[:] / gamma_rot_validate[k, :]))
        drag_steps_0 = int(math.floor(20 * t0_max / self.es.time_step))
        print("drag_steps_0 = {0}".format(drag_steps_0))

        tol = 7E-3
        if "EXTERNAL_FORCES" in espressomd.features():
            for k in range(2):
                self.es.part[k].pos = np.array([0.0, 0.0, 0.0])
                self.es.part[k].v = np.array([0.0, 0.0, 0.0])
                self.es.part[k].omega_body = np.array([0.0, 0.0, 0.0])
            f0 = np.array([-1.2, 58.3578, 0.002])
            f1 = np.array([-15.112, -2.0, 368.0])
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
                    tmp_axis0 = np.cross(tor0, dip0) / (np.linalg.norm(tor0) * np.linalg.norm(dip0))
                    tmp_axis1 = np.cross(tor1, dip1) / (np.linalg.norm(tor1) * np.linalg.norm(dip1))
            self.es.integrator.run(drag_steps_0)
            self.es.time = 0.0
            for k in range(2):
                self.es.part[k].pos = np.array([0.0, 0.0, 0.0])
            if "DIPOLES" in espressomd.features():
                    self.es.part[0].dip = dip0
                    self.es.part[1].dip = dip1
            for i in range(100):
                self.es.integrator.run(10)
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
                    sgn0 = np.sign(np.dot(self.es.part[0].dip, tmp_axis0))
                    sgn0_test = np.sign(np.sin(self.es.time * np.linalg.norm(tor0) / gamma_rot_validate[0, 0]))
                    
                    cos_alpha1 = np.dot(dip1,self.es.part[1].dip) / (np.linalg.norm(dip1) * self.es.part[1].dipm)
                    cos_alpha1_test = np.cos(self.es.time * np.linalg.norm(tor1) / gamma_rot_validate[1, 0])
                    sgn1 = np.sign(np.dot(self.es.part[1].dip, tmp_axis1))
                    sgn1_test = np.sign(np.sin(self.es.time * np.linalg.norm(tor1) / gamma_rot_validate[1, 0]))
                    
                    #print("cos_alpha0 = {0}, cos_alpha0_test={1}".format(cos_alpha0, cos_alpha0_test))
                    #print("dip0 = {0}, self.es.part[0].dip={1}".format(dip0, self.es.part[0].dip))
                    self.assertLess(abs(cos_alpha0 - cos_alpha0_test), tol)
                    self.assertLess(abs(cos_alpha1 - cos_alpha1_test), tol)
                    self.assertEqual(sgn0, sgn0_test)
                    self.assertEqual(sgn1, sgn1_test)

        for i in range(len(self.es.part)):
            self.es.part[i].remove()

        # thermalization
        # Checks if every degree of freedom has 1/2 kT of energy, even when
        # mass and inertia tensor are active

        # 2 different langevin parameters for particles
        temp = np.array([2.5, 2.0])
        D_tr = np.zeros((2, 3))
        D_rot = np.zeros((2, 3))
        for k in range(2):
            gamma_tran[k, :] = np.array((0.4 + np.random.random(3)) * 10)
            gamma_rot[k, :] = np.array((0.2 + np.random.random(3)) * 20)

        box = 10.0
        self.es.box_l = [box, box, box]
        if espressomd.has_features(("PARTIAL_PERIODIC",)):
            self.es.periodicity = 0, 0, 0
        # Random temperature
        kT = (0.3 + np.random.random()) * 5
        gamma_global = np.array((0.5 + np.random.random(3)) * 2.0 / 3.0)
        gamma_global_rot = np.array((0.2 + np.random.random(3)) * 20)

        if test_case in [2,3,(6+self.rot_flag),(7+self.rot_flag)]:
            halfkT = temp / 2.0
        else:
            halfkT = np.array([kT, kT]) / 2.0

        if test_case in [1,3,(5+self.rot_flag),(7+self.rot_flag)]:
            gamma_tr = gamma_tran
        else:
            for k in range(2):
                gamma_tr[k, :] = gamma_global[:]

        if test_case in [1,3,(5+self.rot_flag),(7+self.rot_flag)]:
            gamma_tr = gamma_tran
            gamma_rot_validate = gamma_rot
        else:
            for k in range(2):
                gamma_tr[k, :] = gamma_global[:]
                if (test_case == 4 or test_case == 9) and self.rot_flag == 1:
                    gamma_rot_validate[k, :] = gamma_global_rot[:]
                else:
                    gamma_rot_validate[k, :] = gamma_global[:]

        # translational and rotational diffusion
        for k in range(2):
            D_tr[k, :] = 2.0 * halfkT[k] / gamma_tr[k, :]
            D_rot[k, :] = 2.0 * halfkT[k] / gamma_rot_validate[k, :]

        self.set_thermo_all(test_case, kT, gamma_global, gamma_global_rot)

        # no need to rebuild Verlet lists, avoid it
        self.es.cell_system.skin = 5.0
        if test_case < 4 + self.rot_flag:
            self.es.time_step = 0.03
        else:
            self.es.time_step = 10
        n = 200
        mass = (0.2 + np.random.random()) * 7.0
        J = np.zeros((2, 3))
        for k in range(2):
            J[k, :] = np.array((0.2 + np.random.random(3)) * 7.0)

        for i in range(n):
            for k in range(2):
                ind = i + k * n
                part_pos = np.array(np.random.random(3) * box)
                part_v = np.array([0.0, 0.0, 0.0])
                part_omega_body = np.array([0.0, 0.0, 0.0])
                self.es.part.add(rotation=(1,1,1), id=ind, mass=mass,
                                 rinertia=[J[k, 0], J[k, 1], J[k, 2]],
                                 pos=part_pos, v=part_v)
                if "ROTATION" in espressomd.features():
                    self.es.part[ind].omega_body = part_omega_body
                if test_case in [1,(5+self.rot_flag)]:
                    self.es.part[ind].gamma = gamma_tran[k, :]
                    if "ROTATION" in espressomd.features():
                        self.es.part[ind].gamma_rot = gamma_rot[k, :]

                if test_case in [2,(6+self.rot_flag)]:
                    self.es.part[ind].temp = temp[k]
                if test_case in [3,(7+self.rot_flag)]:
                    self.es.part[ind].gamma = gamma_tran[k, :]
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
        
        # Total curve within a spherical trigonometry.
        # [particle_index, which principal axis, around which lab axis]
        alpha = np.zeros((2, n, 3, 3))
        alpha_norm = np.zeros((2))
        sigma2_alpha = np.zeros((2))
        alpha2 = np.zeros((2))
        # Previous directions of the principal axes:
        # [particle_index, which principal axis, its lab coordinate]
        prev_pa_lab = np.zeros((2, n, 3, 3))
        pa_body = np.zeros((3))
        pa_lab = np.zeros((3, 3))
        ref_lab = np.zeros((3))
        vec = np.zeros((3))
        vec1 = np.zeros((3))
        vec2 = np.zeros((3))
        #vec_diag = np.ones((3))

        pos0 = np.zeros((2 * n, 3))
        for p in range(n):
            for k in range(2):
                ind = p + k * n
                pos0[ind, :] = self.es.part[ind].pos
        dt0 = mass / gamma_tr
        dt0_rot = J / gamma_rot_validate

        #if test_case in [3,(7 + self.rot_flag)]:
        #   loops = 5000
        #else:
        #    loops = 1250
        loops = 5000
        print("Thermalizing...")
        therm_steps = 150
        self.es.integrator.run(therm_steps)
        print("Measuring...")

        int_steps = 1
        fraction_i = 0.65
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
                                                           j] * (2.0 * dt + dt0[k,
                                                                              j] * (- 3.0 + 4.0 * math.exp(- dt / dt0[k,
                                                                                                                  j]) - math.exp(- 2.0 * dt / dt0[k,
                                                                                                                                                j])))
                    dr_norm[k] = dr_norm[k] + \
                        (sum(dr2[k, :]) - sigma2_tr[k]) / sigma2_tr[k]
                    
                    # Rotational diffusion variance.
                    if i >= fraction_i * loops:
                        # let's limit test cases to speed this test..
                        #if test_case in [3,(7 + self.rot_flag)]:
                        dt -= self.es.time_step * (1 + therm_steps + fraction_i * loops)
                        # First, let's identify principal axes in the lab reference frame.
                        alpha2[k] = 0.0
                        sigma2_alpha[k] = 0.0
                        for j in range(3):
                            for j1 in range(3):
                                pa_body[j1] = 0.0
                            pa_body[j] = 1.0
                            vec = tc.convert_vec_body_to_space(self.es, ind, pa_body)
                            pa_lab[j, :] = vec[:]
                            
                            if i >= fraction_i * loops + 1:
                                # Around which axis we rotates?
                                for j1 in range(3):
                                    # Calc a rotational diffusion within the spherical trigonometry
                                    vec2 = vec
                                    vec1[:] = prev_pa_lab[k, p, j, :]
                                    for j2 in range(3):
                                        ref_lab[j2] = 0.0
                                    ref_lab[j1] = 1.0
                                    dalpha = np.arccos(np.dot(vec2, vec1) / (np.linalg.norm(vec2) * np.linalg.norm(vec1)))
                                    rot_projection = np.dot(np.cross(vec2, vec1), ref_lab) / np.linalg.norm(np.cross(vec2, vec1))
                                    theta = np.arccos(np.dot(vec2, ref_lab) / (np.linalg.norm(vec2) * np.linalg.norm(ref_lab)))
                                    alpha[k, p, j, j1] += dalpha * rot_projection / np.sin(theta)
                                    alpha2[k] += alpha[k, p, j, j1]**2
                                    sigma2_alpha[k] += D_rot[k, j] * (2.0 * dt + dt0_rot[k, j] * (- 3.0 +
                                                                                    4.0 * math.exp(- dt / dt0_rot[k, j]) 
                                                                                    - math.exp(- 2.0 * dt / dt0_rot[k, j])))
                            prev_pa_lab[k, p, j, :] = pa_lab[j, :]
                        if i >= fraction_i * loops + 3:
                            alpha_norm[k] += (alpha2[k] - sigma2_alpha[k]) / sigma2_alpha[k]

        tolerance = 0.15
        Ev = 0.5 * mass * v2 / (n * loops)
        Eo = 0.5 * J * o2 / (n * loops)
        dv = np.zeros((2))
        do = np.zeros((2))
        do_vec = np.zeros((2, 3))
        for k in range(2):
            dv[k] = sum(Ev[k, :]) / (3.0 * halfkT[k]) - 1.0
            do[k] = sum(Eo[k, :]) / (3.0 * halfkT[k]) - 1.0
            do_vec[k, :] = Eo[k, :] / halfkT[k] - 1.0
        dr_norm = dr_norm / (n * loops)
        alpha_norm = alpha_norm / (n * (1 - fraction_i) * loops - 2)
        
        for k in range(2):
            print("\n")
            print("k = " + str(k))
            print("mass = " + str(mass))
            print("gamma_tr = {0} {1} {2}".format(
                gamma_tr[k, 0], gamma_tr[k, 1], gamma_tr[k, 2]))
            if test_case in [1,5,3,7]:
                print("gamma_rot_validate = {0} {1} {2}".format(
                    gamma_rot_validate[k, 0], gamma_rot_validate[k, 1], gamma_rot_validate[k, 2]))
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
            print(
                "Deviation in rotational diffusion: {0} ".format(
                    alpha_norm))

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
            if test_case in (1, 3):
                self.assertLessEqual(
                    abs(
                        alpha_norm[k]),
                    tolerance,
                    msg='Relative deviation in rotational diffusion is too large: {0}'.format(
                        alpha_norm[k]))

    def test(self):
        if "ROTATION" in espressomd.features():
            self.rot_flag = 1
        if not "BROWNIAN_DYNAMICS" in espressomd.features():
            n_test_cases = 4 + self.rot_flag
        else:
            n_test_cases = 2 * (4 + self.rot_flag)
        for i in range(n_test_cases):
            self.run_test_case(i)


if __name__ == '__main__':
    print("Features: ", espressomd.features())
    ut.main()
