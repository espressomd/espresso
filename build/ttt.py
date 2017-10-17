#!/usr/bin/env python
 
import espressomd as md
from espressomd import lb
import sys

system = md.System()

system.box_l = [9, 9, 9]
system.time_step = 1.0
system.cell_system.skin = 0.4

lb = md.lb.LBFluid_GPU(agrid=1.0, dens=1.0, visc=1.0, fric=1.0, tau=1.0)
system.actors.add(lb)

while system.time < 1000.1:
  print(lb[0,0,0].velocity)
  sys.stderr.write('t={}\n'.format(system.time))
  system.integrator.run(1)
  system.lees_edwards_offset += 0.01
