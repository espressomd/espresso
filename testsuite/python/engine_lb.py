from __future__ import print_function
import unittest as ut
import tests_common
import os
import numpy as np
import espressomd
from espressomd import lb
vtk_found=True
try:
    import vtk
except:
    vtk_found=False

if "ENGINE" in espressomd.features() and "LB" in espressomd.features():
    class SwimmerTest(ut.TestCase):
        if (vtk_found):
            def test(self):
                # Set to true if you need a new
                # comparison configuration

                new_configuration = False

                boxl      = 12
                sampsteps = 2000
                tstep     = 0.01
                temp      = 0.0
    
                S = espressomd.System()
    
                if (S.cell_system.get_state()["n_nodes"] > 1):
                    print("NOTE: Ignoring testcase for n_nodes > 1")
                    return
    
                S.box_l     = [boxl, boxl, boxl]
                S.cell_system.skin      = 0.1
                S.time_step = tstep
    
                S.part.add(id       = 0, pos=[6.0,3.0,2.0],
                           swimming = {"mode": "pusher", "v_swim": 0.10, "dipole_length": 1.0, "rotational_friction":  2.0},
                           quat     = [np.sqrt(.5),np.sqrt(.5),          0,          0])
                S.part.add(id       = 1, pos=[2.0,3.0,6.0],
                           swimming = {"mode": "pusher", "f_swim": 0.03, "dipole_length": 2.0, "rotational_friction": 20.0},
                           quat     = [np.sqrt(.5),          0,np.sqrt(.5),          0])
                S.part.add(id       = 2, pos=[3.0,2.0,6.0],
                           swimming = {"mode": "puller", "v_swim": 0.15, "dipole_length": 0.5, "rotational_friction": 15.0},
                           quat     = [np.sqrt(.5),          0,          0,np.sqrt(.5)])
                S.part.add(id       = 3, pos=[3.0,6.0,2.0],
                           swimming = {"mode": "puller", "f_swim": 0.05, "dipole_length": 1.5, "rotational_friction":  6.0},
                           quat     = [          0,          0,np.sqrt(.5),np.sqrt(.5)])
    
                lbm = lb.LBFluid(agrid=1.0, tau=tstep, fric=0.5, visc=1.0, dens=1.0)
                S.actors.add(lbm)
    
                #thermostat lb $temp
    
                S.integrator.run(sampsteps)
    
                if new_configuration:
                    lbm.print_vtk_velocity("engine_lb.vtk")
                    self.assertTrue( True )
                else:
                    lbm.print_vtk_velocity("engine_lb_tmp.vtk")
                    different, difference = tests_common.calculate_vtk_max_pointwise_difference("engine_lb.vtk", "engine_lb_tmp.vtk",tol=2.0e-7)
                    os.remove("engine_lb_tmp.vtk")
                    print("Maximum deviation to the reference point is: {}".format(difference))
                    self.assertTrue( different )

if __name__ == '__main__':
    print("Features: ", espressomd.features())
    ut.main()
