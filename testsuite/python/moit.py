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
import scipy.spatial.transform as sst

import espressomd


class TestMomentOfInertia(ut.TestCase):
    """
    Check rotation dynamics with masses.

    """
    TIME_STEP = 0.001
    EXT_TORQUE = np.array([2.0, 2.5, 3.0])
    MOMENTS_OF_INERTIA = np.array([2.4, 3.2, 4.1])

    @classmethod
    def setUpClass(cls):
        cls.system = espressomd.System(box_l=[10.] * 3)
        cls.system.cell_system.skin = 0.4
        cls.system.time_step = cls.TIME_STEP

    def tearDown(self):
        self.system.part.clear()

    @utx.skipIfMissingFeatures(
        ["EXTERNAL_FORCES", "ROTATION", "ROTATIONAL_INERTIA"])
    def test_angular_momentum(self):
        """
        Test that the differential equation for rotation is fulfilled.

        .. math:: \vec{omega}_i = \frac{T_i}{I_i} t

        omega: angular velocity
        T: torque
        I: moment of inertia
        t: time

        """
        part = self.system.part.add(
            pos=0.5 *
            self.system.box_l,
            quat=[
                1,
                0,
                0,
                0],
            ext_torque=self.EXT_TORQUE,
            rotation=[
                1,
                1,
                1],
            rinertia=self.MOMENTS_OF_INERTIA)

        def check(system, part, rotation):
            """
            Verify that angular momentum increases linearly in time.

            Parameter
            ---------
            system : instance of :class:`espressomd.system.System`
                System to integrate.
            part : instance of :class:`espressomd.particle_data.ParticleHandle`
                Particle to check.
            rotation : instance of :class:`scipy.spatial.transform.Rotation`
                Define the particle orientation by a scipy rotation.

            """
            part.omega_body = np.zeros(3)
            # np.roll needed because scipy has real part last.
            part.quat = np.roll(rotation.as_quat(), 1)
            for i in range(10):
                system.integrator.run(1)
                expected = (i + 1) * system.time_step * \
                    rotation.inv().apply(part.ext_torque) / np.copy(part.rinertia)
                np.testing.assert_almost_equal(
                    part.omega_body, expected, decimal=6)

        check(self.system, part, sst.Rotation.from_quat([0, 0, 0, 1]))
        check(self.system, part, sst.Rotation.from_quat(
            [1. / np.sqrt(2), 0, 0, 1. / np.sqrt(2)]))
        check(self.system, part, sst.Rotation.from_quat(
            [0, 1. / np.sqrt(2), 0, 1. / np.sqrt(2)]))
        check(self.system, part, sst.Rotation.from_quat(
            [0, 0, 1. / np.sqrt(2), 1. / np.sqrt(2)]))


if __name__ == "__main__":
    ut.main()
