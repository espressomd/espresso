
#
# Copyright (C) 2017 The ESPResSo project
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
# Tests particle property setters/getters
from __future__ import print_function
import unittest as ut
import espressomd
import numpy as np
from espressomd.electrostatics import *
from espressomd import scafacos

@ut.skipIf(not espressomd.has_features(["ELECTROSTATICS"]),
           "Features not available, skipping test!")
class CoulombCloudWallTune(ut.TestCase):
    """This compares p3m, p3m_gpu electrostatic forces against stored data."""
    S = espressomd.System()
    forces = {}
    tolerance = 1E-3

    def setUp(self):
        self.S.box_l = (10, 10, 10)
        self.S.time_step = 0.01
        self.S.cell_system.skin = 0.4

        #  Clear actors that might be left from prev tests
        if len(self.S.actors):
            del self.S.actors[0]
        self.S.part.clear()
        data = np.genfromtxt("data/coulomb_cloud_wall_system.data")

        # Add particles to system and store reference forces in hash
        # Input format: id pos q f
        for particle in data:
            id = particle[0]
            pos = particle[1:4]
            q = particle[4]
            f = particle[5:]
            self.S.part.add(id=int(id), pos=pos, q=q)
            self.forces[id] = f

    def compare(self, method_name):
        # Compare forces now in the system to stored ones
        force_abs_diff = 0.
        for p in self.S.part:
            force_abs_diff += abs(np.sqrt(sum((p.f - self.forces[p.id])**2)))
        force_abs_diff /= len(self.S.part)

        print(method_name, "force difference", force_abs_diff)

        self.assertTrue(force_abs_diff <= self.tolerance, "Asbolute force difference " +
                        str(force_abs_diff) + " too large for method " + method_name)

    # Tests for individual methods
    if espressomd.has_features(["P3M"]):
        def test_p3m(self):
            # We have to add some tolerance here, because the reference
            # system is not homogeneous
            self.S.actors.add(P3M(bjerrum_length=1., accuracy = 5e-4,
                                  tune=True))
            self.S.integrator.run(0)
            self.compare("p3m")

    if espressomd.has_features(["ELECTROSTATICS", "CUDA"]):
        def test_p3m_gpu(self):
            # We have to add some tolerance here, because the reference
            # system is not homogeneous
            self.S.actors.add(P3M_GPU(bjerrum_length=1., accuracy = 5e-4,
                                      tune=True))
            self.S.integrator.run(0)
            self.compare("p3m_gpu")


if __name__ == "__main__":
    ut.main()
