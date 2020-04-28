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
import tests_common


@utx.skipIfMissingFeatures(["MASS", "ROTATIONAL_INERTIA"])
class RotationalInertia(ut.TestCase):
    longMessage = True
    # Handle for espresso system
    system = espressomd.System(box_l=[1.0, 1.0, 1.0])
    system.cell_system.skin = 0
    # Particle's angular momentum: initial and ongoing
    L_0_lab = np.zeros((3))
    L_lab = np.zeros((3))

    # Angular momentum
    def L_body(self, part):
        return self.system.part[part].omega_body[:] * \
            self.system.part[part].rinertia[:]

    # Set the angular momentum
    def set_L_0(self, part):
        L_0_body = self.L_body(part)
        self.L_0_lab = tests_common.convert_vec_body_to_space(
            self.system, part, L_0_body)

    def set_L(self, part):
        L_body = self.L_body(part)
        self.L_lab = tests_common.convert_vec_body_to_space(
            self.system, part, L_body)

    def test_stability(self):
        self.system.part.clear()
        self.system.part.add(
            pos=np.array([0.0, 0.0, 0.0]), id=0, rotation=(1, 1, 1))

        # Inertial motion around the stable and unstable axes

        tol = 4E-3
        # Anisotropic inertial moment. Stable axes correspond to J[1] and J[2].
        # The unstable axis corresponds to J[0]. These values relation is J[1]
        # < J[0] < J[2].
        J = np.array([5, 0.5, 18.5])

        self.system.part[0].rinertia = J[:]
        # Validation of J[1] stability
        # ----------------------------
        self.system.time_step = 0.0006
        # Stable omega component should be larger than other components.
        stable_omega = 57.65
        self.system.part[0].omega_body = np.array([0.15, stable_omega, -0.043])
        self.set_L_0(0)

        for i in range(100):
            self.set_L(0)
            for k in range(3):
                self.assertAlmostEqual(
                    self.L_lab[k], self.L_0_lab[k], delta=tol,
                    msg='Inertial motion around stable axis J1: Deviation in '
                        'angular momentum is too large. Step {0}, coordinate '
                        '{1}, expected {2}, got {3}'.format(
                        i, k, self.L_0_lab[k], self.L_lab[k]))
            self.assertAlmostEqual(
                self.system.part[0].omega_body[1], stable_omega, delta=tol,
                msg='Inertial motion around stable axis J1: Deviation in omega '
                    'is too large. Step {0}, coordinate 1, expected {1}, got {2}'
                    .format(i, stable_omega, self.system.part[0].omega_body[1]))
            self.system.integrator.run(10)

        # Validation of J[2] stability
        # ----------------------------
        self.system.time_step = 0.01
        # Stable omega component should be larger than other components.
        stable_omega = 3.2
        self.system.part[0].omega_body = np.array(
            [0.011, -0.043, stable_omega])
        self.set_L_0(0)

        for i in range(100):
            self.set_L(0)
            for k in range(3):
                self.assertAlmostEqual(
                    self.L_lab[k], self.L_0_lab[k], delta=tol,
                    msg='Inertial motion around stable axis J2: Deviation in '
                        'angular momentum is too large. Step {0}, coordinate '
                        '{1}, expected {2}, got {3}'.format(
                        i, k, self.L_0_lab[k], self.L_lab[k]))
            self.assertAlmostEqual(
                self.system.part[0].omega_body[2], stable_omega, delta=tol,
                msg='Inertial motion around stable axis J2: Deviation in omega '
                    'is too large. Step {0}, coordinate 2, expected {1}, got {2}'
                    .format(i, stable_omega, self.system.part[0].omega_body[2]))
            self.system.integrator.run(10)

        # Validation of J[0]
        # ------------------
        self.system.time_step = 0.001
        # Unstable omega component should be larger than other components.
        unstable_omega = 5.76
        self.system.part[0].omega_body = np.array(
            [unstable_omega, -0.043, 0.15])
        self.set_L_0(0)

        for i in range(100):
            self.set_L(0)
            for k in range(3):
                self.assertAlmostEqual(
                    self.L_lab[k], self.L_0_lab[k], delta=tol,
                    msg='Inertial motion around stable axis J0: Deviation in '
                        'angular momentum is too large. Step {0}, coordinate '
                        '{1}, expected {2}, got {3}'.format(
                        i, k, self.L_0_lab[k], self.L_lab[k]))
            self.system.integrator.run(10)

    def energy(self, p):
        return 0.5 * np.dot(p.rinertia, p.omega_body**2)

    def momentum(self, p):
        return np.linalg.norm(p.rinertia * p.omega_body)

    def test_energy_and_momentum_conservation(self):
        system = self.system
        system.part.clear()
        system.thermostat.turn_off()
        p = system.part.add(pos=(0, 0, 0), rinertia=(1.1, 1.3, 1.5),
                            rotation=(1, 1, 1), omega_body=(2, 1, 4))
        E0 = self.energy(p)
        m0 = self.momentum(p)
        system.time_step = 0.001
        for _ in range(1000):
            system.integrator.run(100)
            self.assertAlmostEqual(self.energy(p), E0, places=3)
            self.assertAlmostEqual(self.momentum(p), m0, places=3)


if __name__ == '__main__':
    ut.main()
