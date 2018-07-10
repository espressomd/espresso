from __future__ import print_function

import espressomd as md
from espressomd import lb 

import numpy as np
import unittest as ut 

box_l = 9
sys = md.System(box_l=[box_l, box_l, box_l])
sys.time_step = 0.01
sys.cell_system.skin = 0.4

#@ut.skipIf(not md.has_features(["LB"]) or not md.has_features(["LB_GPU"]),
#          "Features not available, skipping test!")
class LBTest(ut.TestCase):

  def test(self):

    lb = md.lb.LBFluidGPU(agrid=1.0, dens=1.0, visc=1.0, fric=1.0, tau=0.01)
    sys.actors.add(lb)
    sys.integrator.run(100)
    
    bulk_stress_observable = np.copy(lb.pi_bulk)
    bulk_stress_manual = 0.0
    
    for i in range(box_l):
      for j in range(box_l):
        for k in range(box_l):
          bulk_stress_manual += lb[i, j, k].pi
    bulk_stress_manual /= (box_l**3)
    
    np.testing.assert_allclose(bulk_stress_observable, bulk_stress_manual, err_msg='Stresses after initializing are not equal')

    sys.thermostat.set_lb(1.0)
    sys.integrator.run(100)
    
    bulk_stress_observable = np.copy(lb.pi_bulk)
    bulk_stress_manual = 0.0
    
    for i in range(box_l):
      for j in range(box_l):
        for k in range(box_l):
          bulk_stress_manual += (lb[i, j, k].pi)
    bulk_stress_manual /= (box_l**3)
    
    np.testing.assert_allclose(bulk_stress_observable, bulk_stress_manual, err_msg= 'Stresses after integration with thermalized LB are not equal')
    
if __name__ == "__main__":
  ut.main()
