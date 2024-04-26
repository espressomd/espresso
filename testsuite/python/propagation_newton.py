#
# Copyright (C) 2013-2022 The ESPResSo project
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
import espressomd.propagation
import itertools
import numpy as np
import unittest as ut
import unittest_decorators as utx


class VelocityVerlet(ut.TestCase):

    system = espressomd.System(box_l=[10., 10., 10.])
    system.cell_system.skin = 0.4
    system.time_step = 0.01

    def tearDown(self):
        self.system.part.clear()
        self.system.thermostat.turn_off()
        self.system.time_step = 0.01

    def calc_pos(self, p, x0, v0):
        t = self.system.time
        return np.copy(0.5 * p.ext_force / p.mass * t**2 + v0 * t + x0)

    def calc_vel(self, p, v0):
        t = self.system.time
        return np.copy(p.ext_force / p.mass * t + v0)

    def calc_rot(self, p, o0):
        t = self.system.time
        return np.copy(p.ext_torque / p.rinertia * t + o0)

    def test_newton_laws(self):
        """
        Tests integration of Newton's equations for a particle with and without
        external forces and with time step changes on the way.

        """
        system = self.system
        system.integrator.set_vv()

        # Newton's 1st law with time step change on the way
        p = system.part.add(pos=(0, 0, 0), v=(1, 2, 3))
        system.time_step = 0.01
        system.time = 12.
        np.testing.assert_allclose(np.copy(p.v), (1, 2, 3))
        for i in range(10):
            np.testing.assert_allclose(np.copy(p.pos), np.copy(
                i * system.time_step * p.v), atol=1E-12)
            system.integrator.run(1)

        # Check that the time has passed
        np.testing.assert_allclose(system.time, 12. + 10 * system.time_step)

        v0 = np.copy(p.v)
        pos0 = np.copy(p.pos)
        system.time_step = 0.02
        for i in range(10):
            np.testing.assert_allclose(np.copy(p.pos), np.copy(
                pos0 + i * system.time_step * p.v), atol=1E-12)
            system.integrator.run(1)

        # Newton's 2nd law
        if espressomd.has_features("EXTERNAL_FORCES"):
            system.time = 0.
            if espressomd.has_features("MASS"):
                p.mass = 2.3
            p.pos = (0, 0, 0)
            pos0 = np.copy(p.pos)
            ext_force = np.array([-2., 1.3, 1.])
            p.ext_force = ext_force
            system.time_step = 0.03
            for i in range(10):
                ref_pos = self.calc_pos(p, pos0, v0)
                np.testing.assert_allclose(np.copy(p.pos), ref_pos, atol=1E-12)
                system.integrator.run(1)

    @utx.skipIfMissingFeatures("VIRTUAL_SITES_RELATIVE")
    def test_07__virtual(self):
        system = self.system
        system.time_step = 0.01

        virtual = system.part.add(pos=[0., 0., 0.], v=[-1., 0., 0.])
        physical = system.part.add(pos=[0., 0., 0.], v=[1., 0., 0.])
        virtual.vs_relative = (physical.id, 0.3, (1., 0., 0., 0.))

        system.integrator.set_vv()
        system.integrator.run(1)

        np.testing.assert_allclose(np.copy(physical.pos), [0.01, 0., 0.])
        np.testing.assert_allclose(np.copy(virtual.pos), [0.01, 0., 0.3])
        np.testing.assert_allclose(np.copy(physical.v), [1., 0., 0.])
        np.testing.assert_allclose(np.copy(virtual.v), [1., 0., 0.])

    @utx.skipIfMissingFeatures(["MASS",
                                "ROTATIONAL_INERTIA",
                                "EXTERNAL_FORCES"])
    def test_propagation(self):
        """
        Check integration of Newton's equations of motion and Euler's equations
        of rotation for various combinations of propagation modes.
        """
        Propagation = espressomd.propagation.Propagation
        system = self.system
        system.time_step = 0.01
        system.integrator.set_vv()
        pos0 = [0., 1., 2.]
        v0 = np.array([1., 2., 3.])
        o0 = np.array([2., 3., 4.])
        ext_force = np.array([-1., +2., -4.])
        ext_torque = np.array([+1., -3., +5.])
        modes_trans = [
            Propagation.NONE,
            Propagation.SYSTEM_DEFAULT,
            Propagation.TRANS_NEWTON,
        ]
        modes_rot = [
            Propagation.NONE,
            Propagation.SYSTEM_DEFAULT,
            Propagation.ROT_EULER,
        ]
        for mode_trans, mode_rot in itertools.product(modes_trans, modes_rot):
            try:
                rotation = (mode_rot != Propagation.NONE) or \
                           (mode_trans == Propagation.SYSTEM_DEFAULT)
                system.part.add(pos=pos0, v=v0, omega_lab=o0,
                                ext_force=ext_force, ext_torque=ext_torque,
                                rotation=3 * [rotation],
                                mass=0.5, rinertia=[0.5, 0.5, 0.5],
                                propagation=mode_trans | mode_rot)
            except BaseException:
                continue
        assert len(system.part) > 3
        system.time = 0.
        for _ in range(10):
            system.integrator.run(1)
            for p in system.part.all():
                if (p.propagation & (Propagation.SYSTEM_DEFAULT |
                                     Propagation.TRANS_NEWTON)):
                    ref_pos = self.calc_pos(p, pos0, v0)
                    ref_vel = self.calc_vel(p, v0)
                else:
                    ref_pos = pos0
                    ref_vel = v0
                if (p.propagation & (Propagation.SYSTEM_DEFAULT |
                                     Propagation.ROT_EULER)):
                    ref_rot = self.calc_rot(p, o0)
                else:
                    ref_rot = o0
                pos = np.copy(p.pos)
                vel = np.copy(p.v)
                rot = np.copy(p.omega_lab)
                np.testing.assert_allclose(pos, ref_pos, rtol=1e-9)
                np.testing.assert_allclose(vel, ref_vel, rtol=1e-9)
                np.testing.assert_allclose(rot, ref_rot, rtol=1e-9)


if __name__ == "__main__":
    ut.main()
