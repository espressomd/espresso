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
sys.set_random_state_PRNG()

def FluidStressNodes(lb, box_l, fluid_stress_nodes):
    
  for i in range(box_l):
    for j in range(box_l):
      for k in range(box_l):
        fluid_stress_nodes += lb[i, j, k].pi
  fluid_stress_nodes /= (box_l**3)
  
  return fluid_stress_nodes

@ut.skipIf(not md.has_features(["LB"]),
          "Features not available, skipping test!")
class LBTest(ut.TestCase):
  def test(self):

    lb1 = md.lb.LBFluid(agrid=1.0, dens=1.0, visc=1.0, fric=1.0, tau=0.01)
    sys.actors.add(lb1)
    sys.integrator.run(100)
    
    fluid_stress = np.copy(lb1.pi_fluid)
    fluid_stress_nodes = np.zeros((3,3))
    FluidStressNodes(lb1, box_l, fluid_stress_nodes)
    
    np.testing.assert_allclose(fluid_stress, fluid_stress_nodes, err_msg='Stresses after initializing are not equal for the GPU implementation')

    sys.thermostat.set_lb(1.0)
    sys.integrator.run(100)
    
    fluid_stress = np.copy(lb1.pi_fluid)
    fluid_stress_nodes = np.zeros((3,3))
    FluidStressNodes(lb1, box_l, fluid_stress_nodes)

    np.testing.assert_allclose(fluid_stress, fluid_stress_nodes, err_msg= 'Stresses after integration with thermalized LB are not equal for the GPU implementation')
    
    obs = observables.LBFluidStress()
    a = obs.calculate()
    fluid_stress_obs = np.array([[a[0],a[1],a[3]], [a[1],a[2],a[4]], [a[3], a[4], a[5]]])

    np.testing.assert_allclose(fluid_stress_obs, fluid_stress_nodes, err_msg= 'Stresses after integration with thermalized LB are not equal for the GPU implementation using the observable')
    
    sys.actors.remove(lb1)

@ut.skipIf(not md.has_features(["LB_GPU"]) or md.has_features(["SHANCHEN"]) ,
          "Features not available, skipping test!")
class LBTestGPU(ut.TestCase):
  def test(self):

    lb2 = md.lb.LBFluidGPU(agrid=1.0, dens=1.0, visc=1.0, fric=1.0, tau=0.01)
    sys.actors.add(lb2)
    sys.integrator.run(100)
    
    fluid_stress = np.copy(lb2.pi_fluid)
    fluid_stress_nodes = np.zeros((3,3))
    FluidStressNodes(lb2, box_l, fluid_stress_nodes)
    
    np.testing.assert_allclose(fluid_stress, fluid_stress_nodes, err_msg='Stresses after initializing are not equal for the CPU implementation')

    sys.thermostat.set_lb(1.0)
    sys.integrator.run(100)
    
    fluid_stress = np.copy(lb2.pi_fluid)
    fluid_stress_nodes = np.zeros((3,3)) 
    FluidStressNodes(lb2, box_l, fluid_stress_nodes)

    np.testing.assert_allclose(fluid_stress, fluid_stress_nodes, err_msg= 'Stresses after integration with thermalized LB are not equal for the CPU implementation')

    obs = observables.LBFluidStress()
    a = obs.calculate()
    fluid_stress_obs = np.array([[a[0],a[1],a[3]], [a[1],a[2],a[4]], [a[3], a[4], a[5]]])

    np.testing.assert_allclose(fluid_stress_obs, fluid_stress_nodes, err_msg= 'Stresses after integration with thermalized LB are not equal for the CPU implementation using the observable')

    sys.actors.remove(lb2)
      
if __name__ == "__main__":
  ut.main()
