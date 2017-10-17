

#
# Copyright (C) 2013,2014,2015,2016 The ESPResSo project
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
from espressomd.interactions import FeneBond
from tests_common import verify_lj_forces
from numpy import random


@ut.skipIf(not espressomd.has_features("ROTATION"),
           "Test requires ROTATION")
class Rotation(ut.TestCase):
    s = espressomd.System()
    s.cell_system.skin=0
    s.time_step=0.01

    def test_langevin(self):
      """Applies langevin thermostat and checks that correct axes get thermalized"""
      s=self.s
      s.thermostat.set_langevin(gamma=1,kT=1)
      for x in 0,1:
          for y in 0,1:
              for z in 0,1:
                  s.part.clear()
                  s.part.add(id=0,pos=(0,0,0),rotation=(x,y,z),quat=(1,0,0,0),omega_body=(0,0,0),torque_lab=(0,0,0))
                  s.integrator.run(500)
                  self.validate(x,0)
                  self.validate(y,1)
                  self.validate(z,2)
    
    def validate(self,rotate,coord):
        if rotate:
            #self.assertNotEqual(self.s.part[0].torque_body[coord],0)
            self.assertNotEqual(self.s.part[0].omega_body[coord],0)
        else:
            #self.assertEqual(self.s.part[0].torque_body[coord],0)
            self.assertEqual(self.s.part[0].omega_body[coord],0)

    def test_axes_changes(self):
        """Verifies that rotation axes in body and space frame stay the same and other axes dont"""
        s=self.s
        s.part.clear()
        s.part.add(id=0,pos=(0.9,0.9,0.9),ext_torque=(1,1,1))
        s.thermostat.turn_off()
        for dir in 0,1,2:
            # Reset orientation
            s.part[0].quat=1,0,0,0
            
            # Enable rotation in a single direction
            rot=[0,0,0]
            rot[dir]=1
            s.part[0].rotation=rot
            
            s.integrator.run(30)

            s.integrator.run(100)

            # Check other axes:
            for axis in [1,0,0],[0,1,0],[0,0,1]:
                if rot==axis:
                    # The axis for which rotation is on should coincide in body and space frame
                    self.assertAlmostEqual(np.dot(rot,s.part[0].convert_vector_body_to_space(rot)),1,places=8)
                else:
                    # For non-rotation axis, body and space frame should differ
                    self.assertLess(np.dot(axis,s.part[0].convert_vector_body_to_space(axis)),0.95)

                    
                 

        




if __name__ == "__main__":
    ut.main()
