from __future__ import print_function
import unittest as ut
import numpy as np
import espressomd
import math


@ut.skipIf(not espressomd.has_features(["MASS", "ROTATIONAL_INERTIA"]),
           "Features not available, skipping test!")
class RotationalInertia(ut.TestCase):
    longMessage = True
    # Handle for espresso system
    es = espressomd.System(box_l=[1.0, 1.0, 1.0])
    es.cell_system.skin = 0
    es.seed = es.cell_system.get_state()['n_nodes'] * [1234]

    def define_rotation_matrix(self, part):
        A = np.zeros((3, 3))
        quat = self.es.part[part].quat
        qq = np.power(quat, 2)

        A[0, 0] = qq[0] + qq[1] - qq[2] - qq[3]
        A[1, 1] = qq[0] - qq[1] + qq[2] - qq[3]
        A[2, 2] = qq[0] - qq[1] - qq[2] + qq[3]

        A[0, 1] = 2 * (quat[1] * quat[2] + quat[0] * quat[3])
        A[0, 2] = 2 * (quat[1] * quat[3] - quat[0] * quat[2])
        A[1, 0] = 2 * (quat[1] * quat[2] - quat[0] * quat[3])

        A[1, 2] = 2 * (quat[2] * quat[3] + quat[0] * quat[1])
        A[2, 0] = 2 * (quat[1] * quat[3] + quat[0] * quat[2])
        A[2, 1] = 2 * (quat[2] * quat[3] - quat[0] * quat[1])

        return A

    def convert_vec_body_to_space(self, part, vec):
        A = self.define_rotation_matrix(part)
        return np.dot(A.transpose(), vec)

    # Angular momentum
    def L_body(self, part):
        return self.es.part[part].omega_body[:] * \
            self.es.part[part].rinertia[:]

    def test_stability(self):
        self.es.part.clear()
        self.es.part.add(pos=np.array([0.0, 0.0, 0.0]), id=0,rotation=(1,1,1))

        # Inertial motion around the stable and unstable axes

        tol = 4E-3
        # Anisotropic inertial moment. Stable axes correspond to J[1] and J[2].
        # The unstable axis corresponds to J[0]. These values relation is J[1]
        # < J[0] < J[2].
        J = np.array([5, 0.5, 18.5])

        # Validation of J[1] stability
        # ----------------------------
        self.es.time_step = 0.0006
        # Stable omega component should be larger than other components.
        stable_omega = 57.65
        self.es.part[0].omega_body = np.array([0.15, stable_omega, -0.043])
        self.es.part[0].rinertia = J[:]

        # Angular momentum
        L_0_body = self.L_body(0)
        L_0_lab = self.convert_vec_body_to_space(0, L_0_body)

        for i in range(100):
            L_body = self.L_body(0)
            L_lab = self.convert_vec_body_to_space(0, L_body)
            for k in range(3):
                self.assertLessEqual(
                    abs(
                        L_lab[k] -
                        L_0_lab[k]),
                    tol,
                    msg='Inertial motion around stable axis J1: Deviation in angular momentum is too large. Step {0}, coordinate {1}, expected {2}, got {3}'.format(
                        i,
                        k,
                        L_0_lab[k],
                        L_lab[k]))
            self.assertLessEqual(
                abs(
                    self.es.part[0].omega_body[1] -
                    stable_omega),
                tol,
                msg='Inertial motion around stable axis J1: Deviation in omega is too large. Step {0}, coordinate 1, expected {1}, got {2}'.format(
                    i,
                    stable_omega,
                    self.es.part[0].omega_body[1]))
            self.es.integrator.run(10)

        # Validation of J[2] stability
        # ----------------------------
        self.es.time_step = 0.01
        # Stable omega component should be larger than other components.
        stable_omega = 3.2
        self.es.part[0].omega_body = np.array([0.011, -0.043, stable_omega])
        self.es.part[0].rinertia = J[:]

        L_0_body = self.L_body(0)
        L_0_lab = self.convert_vec_body_to_space(0, L_0_body)

        for i in range(100):
            L_body = self.L_body(0)
            L_lab = self.convert_vec_body_to_space(0, L_body)
            for k in range(3):
                self.assertLessEqual(
                    abs(
                        L_lab[k] -
                        L_0_lab[k]),
                    tol,
                    msg='Inertial motion around stable axis J2: Deviation in angular momentum is too large. Step {0}, coordinate {1}, expected {2}, got {3}'.format(
                        i,
                        k,
                        L_0_lab[k],
                        L_lab[k]))
            self.assertLessEqual(
                abs(
                    self.es.part[0].omega_body[2] -
                    stable_omega),
                tol,
                msg='Inertial motion around stable axis J2: Deviation in omega is too large. Step {0}, coordinate 2, expected {1}, got {2}'.format(
                    i,
                    stable_omega,
                    self.es.part[0].omega_body[2]))
            self.es.integrator.run(10)

        # Validation of J[0]
        # ------------------
        self.es.time_step = 0.001
        # Unstable omega component should be larger than other components.
        unstable_omega = 5.76
        self.es.part[0].omega_body = np.array([unstable_omega, -0.043, 0.15])
        self.es.part[0].rinertia = J[:]

        L_0_body = self.L_body(0)
        L_0_lab = self.convert_vec_body_to_space(0, L_0_body)

        for i in range(100):
            L_body = self.L_body(0)
            L_lab = self.convert_vec_body_to_space(0, L_body)
            for k in range(3):
                self.assertLessEqual(
                    abs(
                        L_lab[k] -
                        L_0_lab[k]),
                    tol,
                    msg='Inertial motion around stable axis J0: Deviation in angular momentum is too large. Step {0}, coordinate {1}, expected {2}, got {3}'.format(
                        i,
                        k,
                        L_0_lab[k],
                        L_lab[k]))
            self.es.integrator.run(10)


    
    def energy(self,p):
        return 0.5*np.dot(p.rinertia , p.omega_body**2)
    def momentum(self,p):
        return np.linalg.norm(p.rinertia * p.omega_body)

    def test_energy_and_momentum_conservation(self):
        es=self.es
        es.part.clear()
        es.thermostat.turn_off()
        p=es.part.add(pos=(0,0,0),rinertia=(1.1,1.3,1.5),rotation=(1,1,1),omega_body=(2,1,4))
        E0=self.energy(p)
        m0=self.momentum(p)
        es.time_step=0.001
        for i in range(1000):
            es.integrator.run(100)
            self.assertAlmostEqual(self.energy(p),E0,places=3)
            self.assertAlmostEqual(self.momentum(p),m0,places=3)

        



if __name__ == '__main__':
    print("Features: ", espressomd.features())
    ut.main()
