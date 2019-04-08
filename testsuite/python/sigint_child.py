from __future__ import print_function
import numpy as np
import espressomd
import time

system = espressomd.System(box_l=[100, 100, 100])
system.time_step = 0.01
system.cell_system.skin = 0.1

for i in range(100):
    system.part.add(pos=np.random.random() * system.box_l)

while True:
    system.integrator.run(1000)
