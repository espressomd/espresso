#
# Copyright (C) 2013-2019 The ESPResSo project
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
import unittest as ut
import unittest_decorators as utx
import espressomd


@utx.skipIfMissingFeatures(["NPT", "LENNARD_JONES"])
class IntegratorNPT(ut.TestCase):

    """This tests the NpT integrator interface."""
    system = espressomd.System(box_l=[1.0, 1.0, 1.0])

    def setUp(self):
        self.system.box_l = [5] * 3
        self.system.time_step = 0.01
        self.system.cell_system.skin = 0.25

    def tearDown(self):
        self.system.part.clear()
        self.system.thermostat.turn_off()
        self.system.integrator.set_vv()

    def test_integrator_exceptions(self):
        # invalid parameters should throw exceptions
        with self.assertRaises(RuntimeError):
            self.system.integrator.set_isotropic_npt(ext_pressure=-1, piston=1)
        with self.assertRaises(RuntimeError):
            self.system.integrator.set_isotropic_npt(ext_pressure=1, piston=0)
        with self.assertRaises(RuntimeError):
            self.system.integrator.set_isotropic_npt(ext_pressure=1, piston=-1)
        with self.assertRaises(RuntimeError):
            self.system.integrator.set_isotropic_npt(
                ext_pressure=1, piston=1, direction=[0, 0, 0])
        with self.assertRaises(Exception):
            self.system.integrator.set_isotropic_npt(
                ext_pressure=1, piston=1, direction=[1, 0])


if __name__ == "__main__":
    ut.main()
