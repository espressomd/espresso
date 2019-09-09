# Copyright (C) 2010-2018 The ESPResSo project
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
import numpy as np
import espressomd


class AnalyzeMassRelated(ut.TestCase):

    """Test analysis routines that involve particle mass. E.g., center of mass, intertia tensor, ...
    Checks that virtual sites (which do not have meaningful mass, are skipped

    """

    box_l = 50.0
    system = espressomd.System(box_l=[box_l, box_l, box_l])
    np.random.seed(seed=123)

    @classmethod
    def setUpClass(cls):
        cls.system.cell_system.skin = 0.4
        cls.system.time_step = 0.01
        cls.system.thermostat.turn_off()
        cls.system.part.add(
            pos=np.random.random((10, 3)) * cls.box_l, type=[0] * 10)
        cls.system.part.add(
            pos=np.random.random((10, 3)) * cls.box_l, type=[1] * 10)
        if espressomd.has_features("VIRTUAL_SITES"):
            cls.system.part.add(
                pos=np.random.random((10, 3)) * cls.box_l, type=[0] * 10, virtual=[True] * 10)
        if espressomd.has_features("MASS"):
            cls.system.part[:].mass = 0.5 + \
                np.random.random(len(cls.system.part))

    def i_tensor(self, ids):
        pslice = self.system.part[ids]

        I = np.zeros((3, 3))

        # Center of mass
        com = np.zeros(3)
        for p in pslice:
            com += p.mass * p.pos
        com /= np.sum(pslice.mass)

        # Eqn from
        # https://en.wikipedia.org/wiki/Moment_of_inertia#Inertia_tensor
        for p in pslice:
            I += p.mass * \
                (np.dot(p.pos - com, p.pos - com) * np.identity(3)
                 - np.outer(p.pos - com, p.pos - com))
        return I

    def test_itensor(self):
        # Particles of type 0
        I0 = self.i_tensor(self.system.part.select(
            lambda p: (not p.virtual) and p.type == 0).id)

        np.testing.assert_allclose(
            I0, self.system.analysis.moment_of_inertia_matrix(p_type=0), atol=1E-9)
        # type=1
        I1 = self.i_tensor(self.system.part.select(type=1).id)
        self.assertTrue(
            np.allclose(I1, self.system.analysis.moment_of_inertia_matrix(p_type=1), atol=1E-9))

    def test_center_of_mass(self):
        no_virtual_type_0 = self.system.part.select(
            lambda p: (not p.virtual) and p.type == 0)
        com = np.zeros(3)
        for p in no_virtual_type_0:
            com += p.pos * p.mass
        com /= np.sum(no_virtual_type_0.mass)

        np.testing.assert_allclose(
            com,
            self.system.analysis.center_of_mass(p_type=0))
    
    def test_angularmomentum(self):
        no_virtual_type_0 = self.system.part.select(
            lambda p: (not p.virtual) and p.type == 0)
        am = np.zeros(3)
        for p in no_virtual_type_0:
            am += p.mass * np.cross(p.pos, p.v)

        np.testing.assert_allclose(
            am,
            self.system.analysis.angular_momentum(p_type=0))
        

if __name__ == "__main__":
    ut.main()
