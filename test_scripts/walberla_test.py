#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from espressomd import lb, shapes
import espressomd
import sys
print("\n######################################")
print("\n#      WaLBeRLa Shizzle Drizzle      #")
print("\n######################################")
print("\n\n")


boxl = 32
system = espressomd.System(box_l=[boxl]*3)

visc = 1.34
dens = 1.
agrid = 1.
time_step = 0.01

system.time_step = time_step
system.periodicity = [1, 1, 1]
skin = system.cell_system.skin = 0.1


WALL_OFFSET = agrid
right_wall = espressomd.shapes.Wall(normal=[1, 0, 0], dist=WALL_OFFSET)
left_wall = espressomd.shapes.Wall(normal=[-1, 0, 0], dist=-(boxl - WALL_OFFSET))
system.integrator.run(0)

lbf = espressomd.lb.LBFluidWalberla(agrid=agrid, dens=dens, visc=visc, tau=time_step,
                                    seed=152)
system.actors.add(lbf)
system.thermostat.set_lb(LB_fluid=lbf, gamma=10, seed=127)
system.integrator.run(0)

equil_rounds = 100
equil_steps = 100


for i in range(equil_rounds):
    print("Running cycle ", i, " of ", equil_rounds)
    system.integrator.run(equil_steps)

sys.exit()
