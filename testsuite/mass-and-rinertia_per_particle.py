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
    es = espressomd.System()

    def run_test_case(self, test_case):
        seed(1)
        # Decelleration
        self.es.time_step = 0.007
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
        gamma_rot[0, :] = np.array((0.5 + random(3)) * 2.0 / 3.0)
        
        gamma_tran[1, 0] = np.array(0.5 + np.random.random())
        gamma_tran[1, 1] = gamma_tran[1, 0]
        gamma_tran[1, 2] = gamma_tran[1, 0]
        gamma_rot[1, 0] = np.array((0.5 + np.random.random()) * 2.0 / 3.0)
        gamma_rot[1, 1] = gamma_rot[1, 0]
        gamma_rot[1, 2] = gamma_rot[1, 0]

        if test_case != 5:
            self.es.thermostat.set_langevin(
                kT=0.0,
                gamma=[
                    gamma_global[0],
                    gamma_global[1],
                    gamma_global[2]])
        else:
            self.es.thermostat.set_langevin(
                kT=0.0,
                gamma=[
                    gamma_global[0],
                    gamma_global[1],
                    gamma_global[2]],
                gamma_rotation=[
                    gamma_global_rot[0],
                    gamma_global_rot[1],
                    gamma_global_rot[2]])

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
            
            if test_case == 5:
                if (k == 0):
                    print("------------------------------------------------")
                    print("Test " + str(test_case) + ": no particle specific values")
                    print("------------------------------------------------")
                    # No assignments are needed.

        if test_case == 1 or test_case == 3:
            gamma_tr = gamma_tran
            gamma_rot_validate = gamma_rot
        else:
            for k in range(2):
                gamma_tr[k, :] = gamma_global[:]
                if test_case != 5:
                    gamma_rot_validate[k, :] = gamma_global[:]
                else:
                    gamma_rot_validate[k, :] = gamma_global_rot[:]

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

        for i in range(len(self.es.part)):
            self.es.part[i].remove()

        # thermalization
        # Checks if every degree of freedom has 1/2 kT of energy, even when
        # mass and inertia tensor are active

        # 2 different langevin parameters for particles
        temp = np.array([2.5, 2.0])
        D_tr = np.zeros((2, 3))
        for k in range(2):
            gamma_tran[k, :] = np.array((0.4 + np.random.random(3)) * 10)
            # Second particle should be isotropic for the rotational diffusion test.
            # It is valid within test cases #1 and #3
            if test_case in (1, 3) and k == 1:
                gamma_rot[k, 0] = np.array((0.4 + np.random.random()) * 20)
                gamma_rot[k, 1] = gamma_rot[k, 0]
                gamma_rot[k, 2] = gamma_rot[k, 0]
            else:
                gamma_rot[k, :] = np.array((0.2 + np.random.random(3)) * 20)

        box = 10.0
        self.es.box_l = [box, box, box]
        if espressomd.has_features(("PARTIAL_PERIODIC",)):
            self.es.periodicity = 0, 0, 0
        # Random temperature
        kT = (0.01 + np.random.random()) * 10
        gamma_global = np.array((0.5 + np.random.random(3)) * 2.0 / 3.0)
        gamma_global_rot = np.array((0.5 + np.random.random(3)) * 2.0 / 3.0)

        if test_case == 2 or test_case == 3:
            halfkT = temp / 2.0
        else:
            halfkT = np.array([kT, kT]) / 2.0

        if test_case == 1 or test_case == 3:
            gamma_tr = gamma_tran
        else:
            for k in range(2):
                gamma_tr[k, :] = gamma_global[:]

        # translational diffusion
        for k in range(2):
            D_tr[k, :] = 2.0 * halfkT[k] / gamma_tr[k, :]

        if test_case != 5:
            self.es.thermostat.set_langevin(
                kT=kT,
                gamma=[
                    gamma_global[0],
                    gamma_global[1],
                    gamma_global[2]])
        else:
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

        # no need to rebuild Verlet lists, avoid it
        self.es.cell_system.skin = 5.0
        self.es.time_step = 0.03
        n = 200
        mass = (0.2 + np.random.random()) * 7.0
        J = np.zeros((2, 3))
        J[0, :] = np.array((0.2 + np.random.random(3)) * 7.0)
        # Second particle should be isotropic for the rotational diffusion test.
        # It is valid within test cases #1 and #3
        if test_case in (1, 3):
            J[1, 0] = np.array((0.2 + np.random.random()) * 1.0)
            J[1, 1] = J[1, 0]
            J[1, 2] = J[1, 0]
        else:
            J[1, :] = np.array((0.2 + np.random.random(3)) * 7.0)

        D_rot_p1 = 2.0 * halfkT[1] / gamma_rot[1, 0]

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
                if test_case == 1:
                    self.es.part[ind].gamma = gamma_tran[k, :]
                    if "ROTATION" in espressomd.features():
                        self.es.part[ind].gamma_rot = gamma_rot[k, :]

                if test_case == 2:
                    self.es.part[ind].temp = temp[k]
                if test_case == 3:
                    self.es.part[ind].gamma = gamma_tran[k, :]
                    if "ROTATION" in espressomd.features():
                        self.es.part[ind].gamma_rot = gamma_rot[k, :]
                    self.es.part[ind].temp = temp[k]
                # It is needed for further rotational diffusion validation:
                if k == 1:
                    if i <= n / 3:
                        self.es.part[ind].rotation = 1, 0, 0
                    if i > n / 3 and i <= 2 * n / 3:
                        self.es.part[ind].rotation = 0, 1, 0
                    if i > 2 * n / 3:
                        self.es.part[ind].rotation = 0, 0, 1

        # Get rid of short range calculations if exclusions are on
        # if espressomd.has_features("EXCLUSIONS"):

        # matrices: [2 types of particless] x [3 dimensions X Y Z]
        # velocity^2, omega^2, position^2
        v2 = np.zeros((2, 3))
        o2 = np.zeros((2, 3))
        dr2 = np.zeros((2, 3))
        sigma2_tr = np.zeros((2))
        dr_norm = np.zeros((2))
        
        # Only for the second particle.
        # Total curve within a spherical trigonometry:
        alpha = np.zeros((n, 3))
        alpha_norm = 0.0
        # Previous directions of the principal axes:
        # [particle_index, which principal axis, its lab coordinate]
        prev_pa_lab = np.zeros((n, 3, 3))
        pa_body = np.zeros((3))
        pa_lab = np.zeros((3, 3))
        vec = np.zeros((3))
        vec1 = np.zeros((3))
        vec_diag = np.ones((3))

        pos0 = np.zeros((2 * n, 3))
        for p in range(n):
            for k in range(2):
                ind = p + k * n
                pos0[ind, :] = self.es.part[ind].pos
        dt0 = mass / gamma_tr
        # Corresponding below test will be made only for the second particle:
        dt0_rot_1 = J[1, 0] / gamma_rot[1, 0]

        loops = 1250
        print("Thermalizing...")
        therm_steps = 150
        self.es.integrator.run(therm_steps)
        print("Measuring...")

        int_steps = 1
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
                    
                    # Rotational diffusion variance.
                    if test_case in (1, 3) and k == 1:
                        dt -= self.es.time_step * (1 + therm_steps)
                        # First, let's identify principal axes in the lab reference frame.
                        for j in range(3):
                            for j1 in range(3):
                                pa_body[j1] = 0.0
                            pa_body[j] = 1.0
                            vec = tc.convert_vec_body_to_space(self.es, ind, pa_body)
                            pa_lab[j, :] = vec[:]
                            
                            if i > 0:
                                # Calc a rotational diffusion within the spherical trigonometry
                                vec1[:] = prev_pa_lab[p, j, :]
                                dalpha = np.arccos(np.dot(vec, vec1) / (np.linalg.norm(vec) * np.linalg.norm(vec1)))
                                # just a formal sign keep to distinguish opposite rotations
                                sign = np.sign(np.dot(np.cross(vec, vec1), vec_diag))
                                alpha[p, j] += sign * dalpha
                                alpha2 = alpha[p, j]**2
                                sigma2_alpha = D_rot_p1 * (2.0 * dt + dt0_rot_1 * (- 3.0 +
                                                                                4.0 * math.exp(- dt / dt0_rot_1) 
                                                                                - math.exp(- 2.0 * dt / dt0_rot_1)))
                                if dalpha > 0:
                                    alpha_norm += (alpha2 - sigma2_alpha) / sigma2_alpha
                            prev_pa_lab[p, j, :] = pa_lab[j, :]

        tolerance = 0.15
        Ev = 0.5 * mass * v2 / (n * loops)
        Eo = 0.5 * J * o2 / (n * loops)
        dv = np.zeros((2))
        do = np.zeros((2))
        do_vec = np.zeros((2, 3))
        for k in range(2):
            dv[k] = sum(Ev[k, :]) / (3 * halfkT[k]) - 1.0
            if k == 0:
                do[k] = sum(Eo[k, :]) / (3 * halfkT[k]) - 1.0
                do_vec[k, :] = Eo[k, :] / halfkT[k] - 1.0
            else:
                # Two rotational axes are fixed for the second particle:
                do[k] = sum(Eo[k, :]) / (1 * halfkT[k]) - 1.0
                do_vec[k, :] = Eo[k, :] / (halfkT[k] / 3.0) - 1.0
        dr_norm = dr_norm / (n * loops)
        
        # Only two body axes move around the non-fixed third one.
        alpha_norm = alpha_norm / (2 * n * loops)
        
        for k in range(2):
            print("\n")
            print("k = " + str(k))
            print("mass = " + str(mass))
            print("gamma_tr = {0} {1} {2}".format(
                gamma_tr[k, 0], gamma_tr[k, 1], gamma_tr[k, 2]))
            if test_case == 1 or test_case == 3:
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
                        alpha_norm),
                    tolerance,
                    msg='Relative deviation in rotational diffusion is too large: {0}'.format(
                        alpha_norm))

    def test(self):
        if "ROTATION" in espressomd.features():
            n_test_cases = 5
        else:
            n_test_cases = 4
        for i in range(n_test_cases):
            self.run_test_case(i)


if __name__ == '__main__':
    print("Features: ", espressomd.features())
    ut.main()
