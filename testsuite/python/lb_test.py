from __future__ import print_function
import espressomd
import espressomd.lb

S = espressomd.System()
S.box_l = [16, 16, 16]
S.cell_system.skin = 0.4
S.time_step = 0.01

print("Setup LB")
lb = espressomd.lb.LBFluid(dens=0.5, agrid=1.0, visc=0.8, tau=0.1, fric=1, ext_force=[1e-1,2e-1,3e-1])

# print "Get LB params", lb.getParams()


print("Add LB Actor")
S.actors.add(lb)
print("Get LB params ", lb.get_params())
lb.set_params(dens=1.0)
print("Get LB params 2nd ", lb.get_params())
#lb.set_params(gpu="yes")


print("Place Particles")
S.part.add(pos=[10,10,25])
S.part.add(pos=[10,10,75])


print("Integrate")
for i in range(0,10):
   S.integrator.run(100)
   integrate.integrate(100)
   print("P1: " + str(S.part[0].pos))
   print("P2: " + str(S.part[1].pos))
