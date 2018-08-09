from __future__ import print_function
import espressomd.lb
import espressomd.lbboundaries
import espressomd.shapes
import unittest as ut
import numpy as np

@ut.skipIf(not espressomd.has_features(["VIRTUAL_SITES"]),
           "Features not available, skipping test.")
class LBBoundaryThermoVirtualTest(ut.TestCase):
    """Test slip velocity of boundaries.

       In this simple test add wall with a slip verlocity is
       added and checkeckt if the fluid obtains the same velocity.
    """

    system = espressomd.System(box_l=[10.0, 10.0, 10.0])

    system.time_step = 1.0
    system.cell_system.skin = 0.1

    def tearDown(self):
        self.system.part.clear()

        for a in self.system.actors:
            self.system.actors.remove(a)

    def check_virtual(self,fluid_class):
        s = self.system
        lb_fluid = fluid_class(
            agrid=1.0, dens=1.0, visc=1.0, fric=1.0, tau=1.0)
        s.actors.add(lb_fluid)


        virtual = s.part.add(pos=[0,0,0], virtual=True, v=[1,0,0])
        physical = s.part.add(pos=[0,0,0], virtual=False, v=[1,0,0])

        s.thermostat.set_lb(kT=0, act_on_virtual=False)

        s.integrator.run(1)

        np.testing.assert_almost_equal(np.copy(virtual.f), [0,0,0])
        np.testing.assert_almost_equal(np.copy(physical.f), [-1,0,0])

        s.thermostat.set_lb(kT=0, act_on_virtual=True)

        virtual.v = [1,0,0]
        physical.v = [1,0,0]

        
        s.actors.remove(lb_fluid)
        lb_fluid = fluid_class(
            agrid=1.0, dens=1.0, visc=1.0, fric=1.0, tau=1.0)
        s.actors.add(lb_fluid)
        virtual.pos=physical.pos
        virtual.v=1,0,0
        physical.v=1,0,0
        s.integrator.run(1)

        # The forces are not exactly -1. because the fluid is not at
        # rest anymore because of the previous check.
        np.testing.assert_almost_equal(np.copy(physical.f),np.copy(virtual.f))
        np.testing.assert_almost_equal(np.copy(physical.f), [-1,0,0])
        np.testing.assert_almost_equal(np.copy(virtual.f), [-1,0,0])

    @ut.skipIf(not espressomd.has_features(["LB"]),
               "Features not available, skipping test.")
    def test_lb_cpu(self):
        self.check_virtual(espressomd.lb.LBFluid)

    @ut.skipIf(not espressomd.has_features(["LB_GPU"]),
               "Features not available, skipping test.")
    def test_lb_gpu(self):
        self.check_virtual(espressomd.lb.LBFluidGPU)

if __name__ == "__main__":
    ut.main()

