from __future__ import print_function
import sys
import unittest as ut
import numpy as np
import espressomd


class AnalyzeITensor(ut.TestCase):
    """Test the inertia tensor analysis"""
    
    box_l = 50.0
    system = espressomd.System(box_l = [box_l, box_l, box_l])
    np.random.seed(seed=123)

    @classmethod
    def setUpClass(self):
        self.system.cell_system.skin = 0.4
        self.system.time_step = 0.01
        self.system.thermostat.turn_off()
        for i in range(10):
            self.system.part.add(id=i, pos=np.random.random(3) * self.box_l,type=0)
        for i in range(20,30,1):
            self.system.part.add(id=i, pos=np.random.random(3) * self.box_l,type=1)
        if espressomd.has_features("mass"):
            s.part[:].mass=0.5+np.random.random(20)

    def i_tensor(self,ids):
        pslice=self.system.part[ids]

        I=np.zeros((3,3))

        # Center of mass
        com=np.zeros(3)
        for p in pslice:
            com+=p.mass*p.pos
        com/=np.sum(pslice.mass)


        # Eqn from https://en.wikipedia.org/wiki/Moment_of_inertia#Inertia_tensor
        for p in pslice:
            I+=p.mass*(np.dot(p.pos-com,p.pos-com)*np.identity(3) - np.outer(p.pos-com,p.pos-com))
        return I
        
  
        
    
    def test(self):
        # Partilces of type 0
        I0=self.i_tensor(range(0,10,1))
        
        np.testing.assert_allclose(I0, self.system.analysis.moment_of_inertia_matrix(p_type=0),atol=1E-9)
        # type=1
        I1=self.i_tensor(range(20,30,1))
        self.assertTrue(np.allclose(I1, self.system.analysis.moment_of_inertia_matrix(p_type=1),atol=1E-9))

        
        



if __name__ == "__main__":
    ut.main()
