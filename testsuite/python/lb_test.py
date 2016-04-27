import espressomd
from espressomd.lb import LB_FLUID 
from espressomd import integrate

S=espressomd.System()
S.box_l=[16,16,16]
S.skin=0.4
S.time_step=0.01

print "Setup LB"
lb=LB_FLUID(dens=0.5,agrid=1.0,visc=0.8,tau=0.1)

#print "Get LB params", lb.getParams()


print "Add LB Actor"
S.actors.add(lb)
print "Get LB params ", lb.getParams()
lb.setParams(dens=1.0)
print "Get LB params 2nd ", lb.getParams()
lb.setParams(gpu="yes")


#print "Place Particles"
#S.part[0].pos=[10,10,25]
#S.part[1].pos=[10,10,75]    


#print "Integrate"
#for i in range(0,10):
#    integrate.integrate(100)
#    print "P1: " + str(S.part[0].pos)
#    print "P2: " + str(S.part[1].pos)


