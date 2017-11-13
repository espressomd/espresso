from __future__ import print_function
import unittest as ut
import numpy as np
import espressomd
from espressomd import lb

@ut.skipIf(not espressomd.has_features(["ROTATION", "PARTICLE_ANISOTROPY", "LANGEVIN_PER_PARTICLE", "DIPOLES", "EXTERNAL_FORCES", "LB", "EXCLUSIONS"]),
           "Features not available, skipping test!")
class ArrayPropertyTest(ut.TestCase):
    
    system = espressomd.System()
    system.box_l = [12.0,12.0,12.0]
    system.time_step = 0.01
    system.cell_system.skin = 0.01
    system.part.add(pos=[0,0,0])
    lbf = lb.LBFluid(agrid=0.5, dens=1, visc=1, tau=0.01)
    system.actors.add(lbf)

    def locked_operators(self, v):
        with self.assertRaises(UserWarning):
            v[0] = 0
        with self.assertRaises(UserWarning):
            v += [1,1,1]
        with self.assertRaises(UserWarning):
            v -= [1,1,1]
        with self.assertRaises(UserWarning):
            v *= [1,1,1]
        with self.assertRaises(UserWarning):
            v /= [1,1,1]
        with self.assertRaises(UserWarning):
            v //= [1,1,1]
        with self.assertRaises(UserWarning):
            v %= [1,1,1]
        with self.assertRaises(UserWarning):
            v **= [1,1,1]
        with self.assertRaises(UserWarning):
            v <<= [1,1,1]
        with self.assertRaises(UserWarning):
            v >>= [1,1,1]
        with self.assertRaises(UserWarning):
            v &= [1,1,1]
        with self.assertRaises(UserWarning):
            v |= [1,1,1]
        with self.assertRaises(UserWarning):
            v ^= [1,1,1]

    def set_copy(self, v):
        cpy = np.copy(v)
        cpy[0] = 1234
        self.assertTrue(cpy[0] == 1234)

    def test(self):
        
        # Check for exception for various operators
        # Particle
        self.locked_operators(self.system.part[0].pos)
        self.locked_operators(self.system.part[0].v)
        self.locked_operators(self.system.part[0].f)
        self.locked_operators(self.system.part[0].pos_folded)
        self.locked_operators(self.system.part[0].omega_lab)
        self.locked_operators(self.system.part[0].quat)
        self.locked_operators(self.system.part[0].rotation)
        self.locked_operators(self.system.part[0].omega_body)
        self.locked_operators(self.system.part[0].torque_lab)
        self.locked_operators(self.system.part[0].rinertia)
        self.locked_operators(self.system.part[0].director)
        self.locked_operators(self.system.part[0].gamma_rot)
        self.locked_operators(self.system.part[0].gamma)
        self.locked_operators(self.system.part[0].dip)
        self.locked_operators(self.system.part[0].ext_force)
        self.locked_operators(self.system.part[0].fix)
        self.locked_operators(self.system.part[0].ext_torque)
        self.locked_operators(self.system.part[0].exclusions)

        # LB
        self.locked_operators(self.lbf[0,0,0].velocity)
        self.locked_operators(self.lbf[0,0,0].density)
        self.locked_operators(self.lbf[0,0,0].pi)
        self.locked_operators(self.lbf[0,0,0].pi_neq)
        self.locked_operators(self.lbf[0,0,0].population)
        
        # System
        self.locked_operators(self.system.box_l)
        self.locked_operators(self.system.periodicity)
        
        # Check (allowed) setter
        # Particle 
        self.system.part[0].pos         = [2,2,2]    
        self.assertTrue( (self.system.part[0].pos == [2,2,2]).all() )

        self.system.part[0].v           = [2,2,2]
        self.assertTrue( (self.system.part[0].v == [2,2,2]).all() )

        self.system.part[0].f           = [2,2,2]
        self.assertTrue( (self.system.part[0].f == [2,2,2]).all() )

        self.system.part[0].quat        = [0.5,0.5,0.5,0.5]
        self.assertTrue( (self.system.part[0].quat == [0.5,0.5,0.5,0.5]).all() )

        self.system.part[0].torque_lab  = [2,2,2]
        self.assertTrue( (self.system.part[0].torque_lab == [2,2,2]).all() )
        
        self.system.part[0].omega_lab   = [2,2,2]
        self.assertTrue( (self.system.part[0].omega_lab == [2,2,2]).all() )

        self.system.part[0].rotation    = [1,1,1]    
        self.assertTrue( (self.system.part[0].rotation == [1,1,1]).all() )

        self.system.part[0].omega_body  = [2,2,2]
        self.assertTrue( (self.system.part[0].omega_body == [2,2,2]).all() )

        self.system.part[0].rinertia    = [2,2,2]
        self.assertTrue( (self.system.part[0].rinertia == [2,2,2]).all() )

        self.system.part[0].gamma_rot   = [2,2,2]
        self.assertTrue( (self.system.part[0].gamma_rot == [2,2,2]).all() )

        self.system.part[0].gamma       = [2,2,2] 
        self.assertTrue( (self.system.part[0].gamma == [2,2,2]).all() )

        self.system.part[0].dip         = [2,2,2]
        self.assertTrue( (self.system.part[0].dip == [2,2,2]).all() )

        self.system.part[0].ext_force   = [2,2,2]
        self.assertTrue( (self.system.part[0].ext_force == [2,2,2]).all() )

        self.system.part[0].fix         = [1,1,1]
        self.assertTrue( (self.system.part[0].fix == [1,1,1]).all() )

        self.system.part[0].ext_torque  = [2,2,2]
        self.assertTrue( (self.system.part[0].ext_torque == [2,2,2]).all() )
        
        # System        
        self.system.box_l = [2,2,2]
        self.assertTrue( (self.system.box_l == [2,2,2]).all() )
        
        self.system.periodicity = [1,0,0]
        self.assertTrue( (self.system.periodicity == [1,0,0]).all() )

        
        # Check if copy is settable
        # Particle
        self.set_copy(self.system.part[0].pos)
        self.set_copy(self.system.part[0].pos)       
        self.set_copy(self.system.part[0].v)         
        self.set_copy(self.system.part[0].f)         
        self.set_copy(self.system.part[0].pos_folded)
        self.set_copy(self.system.part[0].omega_lab) 
        self.set_copy(self.system.part[0].quat)      
        self.set_copy(self.system.part[0].rotation)  
        self.set_copy(self.system.part[0].omega_body)
        self.set_copy(self.system.part[0].torque_lab)
        self.set_copy(self.system.part[0].rinertia)  
        self.set_copy(self.system.part[0].director)  
        self.set_copy(self.system.part[0].gamma_rot) 
        self.set_copy(self.system.part[0].gamma)     
        self.set_copy(self.system.part[0].dip)       
        self.set_copy(self.system.part[0].ext_force) 
        self.set_copy(self.system.part[0].fix)       
        self.set_copy(self.system.part[0].ext_torque)
        
        # System
        self.set_copy(self.system.box_l)
        self.set_copy(self.system.periodicity)

if __name__ == "__main__":
    ut.main()
