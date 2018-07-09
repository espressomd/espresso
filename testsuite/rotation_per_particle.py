

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
    s = espressomd.System(box_l=[1.0, 1.0, 1.0])
    s.seed = s.cell_system.get_state()['n_nodes'] * [1234]
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

    @ut.skipIf(not espressomd.has_features("EXTERNAL_FORCES"),"Requires EXTERNAL_FORCES")
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

                    


    def test_frame_conversion_and_rotation(self):
        s=self.s
        s.part.clear()
        p=s.part.add(pos=np.random.random(3),rotation=(1,1,1))
        
        # Space and body frame co-incide?
        np.testing.assert_allclose(p.director,p.convert_vector_body_to_space((0,0,1)),atol=1E-10)

        # Random vector should still co-incide
        v=(1.,5.5,17)
        np.testing.assert_allclose(v,p.convert_vector_space_to_body(v),atol=1E-10)
        np.testing.assert_allclose(v,p.convert_vector_body_to_space(v),atol=1E-10)

        # Particle rotation
        
        p.rotate((1,2,0),np.pi/4)
        # Check angle for director
        self.assertAlmostEqual(np.arccos(np.dot(p.director, (0,0,1))),np.pi/4,delta=1E-10)
        # Check other vector
        v=(5,-7,3)
        v_r=p.convert_vector_body_to_space(v)
        self.assertAlmostEqual(np.dot(v,v),np.dot(v_r,v_r),delta=1e-10)
        np.testing.assert_allclose(p.convert_vector_space_to_body(v_r),v,atol=1E-10)

        # Rotation axis should co-incide
        np.testing.assert_allclose((1,2,0),p.convert_vector_body_to_space((1,2,0)))
        
        
        # Check rotation axis with all elements set
        p.rotate(axis=(-5,2,17),angle=1.)
        v=(5,-7,3)
        v_r=p.convert_vector_body_to_space(v)
        self.assertAlmostEqual(np.dot(v,v),np.dot(v_r,v_r),delta=1e-10)
        np.testing.assert_allclose(p.convert_vector_space_to_body(v_r),v,atol=1E-10)










        # 




        




if __name__ == "__main__":
    ut.main()
