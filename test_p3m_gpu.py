import espressomd
from espressomd import electrostatics
from espressomd import *

import random

S=espressomd.System()

S.time_step=0.01
S.skin=0.5
S.box_l=[100, 100, 100]


S.part[0].pos=5.0, 5.0, 5.0
S.part[0].q = 1

S.part[1].pos=0.0, 0.0, 0.0
S.part[1].q = -1

#tcl.eval("inter coulomb 2.0 p3m gpu tune accuracy 1e-2")
#tcl.eval("inter coulomb 2.0 p3m gpu 0.9 24 5 1.0")
#p3m=electrostatics.P3M_GPU(bjerrum_length=2.0, accuracy=1e-2)
p3m=electrostatics.P3M_GPU(bjerrum_length=2.0, accuracy=1e-2, r_cut=0.5, mesh=[24, 24, 24], cao=5)
print "before add"


S.Actors.add(p3m)
print "after add"

for i in range(10):
    integrate(100)
    print S.part


