from __future__ import print_function

import espressomd as md
from espressomd import lb 
from espressomd import observables

import numpy as np
import unittest as ut 

box_l = 9
sys = md.System(box_l=[box_l, box_l, box_l])
sys.time_step = 0.01
sys.cell_system.skin = 0.4

@ut.skipIf(not md.has_features(["LB"]),
          "Features not available, skipping test!")
class LBTest(ut.TestCase):

  def test(self):

    lb1 = md.lb.LBFluid(agrid=1.0, dens=1.0, visc=1.0, fric=1.0, tau=0.01)
    sys.actors.add(lb1)
    sys.integrator.run(100)
    
    bulk_stress_observable = np.copy(lb1.pi_bulk)
    bulk_stress_manual = 0.0
    
    for i in range(box_l):
      for j in range(box_l):
        for k in range(box_l):
          bulk_stress_manual += lb1[i, j, k].pi
    bulk_stress_manual /= (box_l**3)
    
    np.testing.assert_allclose(bulk_stress_observable, bulk_stress_manual, err_msg='Stresses after initializing are not equal for the GPU implementation')

    sys.thermostat.set_lb(1.0)
    sys.integrator.run(100)
    
    bulk_stress_observable = np.copy(lb1.pi_bulk)
    bulk_stress_manual = 0.0
    
    for i in range(box_l):
      for j in range(box_l):
        for k in range(box_l):
          bulk_stress_manual += (lb1[i, j, k].pi)
    bulk_stress_manual /= (box_l**3)
    
    np.testing.assert_allclose(bulk_stress_observable, bulk_stress_manual, err_msg= 'Stresses after integration with thermalized LB are not equal for the GPU implementation')

    obs = observables.LBBulkStress()
    a = obs.calculate()
    bulk_stress_obs = np.array([[a[0],a[1],a[3]], [a[1],a[2],a[4]], [a[3], a[4], a[5]]])

    np.testing.assert_allclose(bulk_stress_obs, bulk_stress_manual, err_msg= 'Stresses after integration with thermalized LB are not equal for the GPU implementation using the observable')
    
    sys.actors.remove(lb1)

@ut.skipIf(not md.has_features(["LB_GPU"]),
          "Features not available, skipping test!")
class LBTestGPU(ut.TestCase):

  def test(self):

    lb2 = md.lb.LBFluidGPU(agrid=1.0, dens=1.0, visc=1.0, fric=1.0, tau=0.01)
    sys.actors.add(lb2)
    sys.integrator.run(100)
    
    bulk_stress_observable = np.copy(lb2.pi_bulk)
    bulk_stress_manual = 0.0
    
    for i in range(box_l):
      for j in range(box_l):
        for k in range(box_l):
          bulk_stress_manual += lb2[i, j, k].pi
    bulk_stress_manual /= (box_l**3)
    
    np.testing.assert_allclose(bulk_stress_observable, bulk_stress_manual, err_msg='Stresses after initializing are not equal for the CPU implementation')

    sys.thermostat.set_lb(1.0)
    sys.integrator.run(100)
    
    bulk_stress_observable = np.copy(lb2.pi_bulk)
    bulk_stress_manual = 0.0
    
    for i in range(box_l):
      for j in range(box_l):
        for k in range(box_l):
          bulk_stress_manual += (lb2[i, j, k].pi)
    bulk_stress_manual /= (box_l**3)
    
    np.testing.assert_allclose(bulk_stress_observable, bulk_stress_manual, err_msg= 'Stresses after integration with thermalized LB are not equal for the CPU implementation')

    obs = observables.LBBulkStress()
    a = obs.calculate()
    bulk_stress_obs = np.array([[a[0],a[1],a[3]], [a[1],a[2],a[4]], [a[3], a[4], a[5]]])

    np.testing.assert_allclose(bulk_stress_obs, bulk_stress_manual, err_msg= 'Stresses after integration with thermalized LB are not equal for the CPU implementation using the observable')

    sys.actors.remove(lb2)
      
if __name__ == "__main__":
  ut.main()
