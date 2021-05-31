#
# Copyright (C) 2020 The ESPResSo project
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
#
import espressomd
import numpy as np
import unittest as ut
import unittest_decorators as utx


@utx.skipIfMissingFeatures("ROTATION")
class Integrate(ut.TestCase):

    system = espressomd.System(box_l=[1., 1., 1.])
    system.cell_system.skin = 0

    def tearDown(self):
        self.system.part.clear()

    @utx.skipIfMissingFeatures("EXTERNAL_FORCES")
    def test_torque_driven_anisotropic(self):
        """
        Apply a constant external torque on 1 axis in the anisotropic case.
        Angular velocity increases linearly.
        """
        system = self.system
        system.time_step = time_step = 0.01
        system.integrator.set_vv()
        # set up 3 particles accelerating around the main axes
        for i in range(3):
            torque = [0, 0, 0]
            torque[i] = 6
            system.part.add(pos=(0, 0, 0), rotation=[True, True, True],
                            ext_torque=torque, rinertia=[1, 2, 3])
        # check angular velocity increases linearly
        for j in range(1, 200):
            system.integrator.run(1)
            for i, p in enumerate(system.part):
                ref_omega = [0, 0, 0]
                ref_omega[i] = j * time_step * p.torque_lab[i] / p.rinertia[i]
                val_omega = np.copy(p.omega_lab)
                np.testing.assert_allclose(val_omega, ref_omega, atol=1e-12)

    @utx.skipIfMissingFeatures("EXTERNAL_FORCES")
    def test_torque_driven_isotropic(self):
        """
        Apply a constant external torque on 3 axes in the isotropic case.
        Angular velocity increases linearly.
        """
        system = self.system
        system.time_step = time_step = 0.01
        system.integrator.set_vv()
        # set up 1 particle accelerating around the 3 main axes
        torque = [6, 7, 8]
        p = system.part.add(pos=(0, 0, 0), rotation=[True, True, True],
                            ext_torque=torque, rinertia=[1, 1, 1])
        # check angular velocity increases linearly
        for j in range(1, 200):
            system.integrator.run(1)
            ref_omega = j * time_step * np.copy(p.torque_lab / p.rinertia)
            val_omega = np.copy(p.omega_lab)
            np.testing.assert_allclose(val_omega, ref_omega, atol=1e-12)

    def test_torque_free_anisotropic(self):
        """
        Rotate particles in the torque-free case.
        """
        system = self.system
        n_steps = 400
        system.time_step = time_step = np.pi / n_steps
        system.integrator.set_vv()
        # set up 3 particles rotating around the main axes
        for i in range(3):
            omega = [0, 0, 0]
            omega[i] = 1
            system.part.add(pos=(0, 0, 0), rotation=[True, True, True],
                            omega_lab=omega, rinertia=[1, 2, 3])
        # check angular velocity remains constant
        # check quaternion follows a sine (the quaternion period is 4 pi)
        for j in range(1, 4 * n_steps + 1):
            system.integrator.run(1)
            for i, p in enumerate(system.part):
                ref_omega = [0, 0, 0]
                ref_omega[i] = 1
                val_omega = np.copy(p.omega_lab)
                np.testing.assert_allclose(val_omega, ref_omega, atol=1e-12)
                ref_quat = np.sin(j * time_step / 2)
                tol = 1e-6 + 1.1e-8 * j  # envelope of the periodic drift
                np.testing.assert_allclose(p.quat[1 + i], ref_quat, atol=tol)


if __name__ == "__main__":
    ut.main()
