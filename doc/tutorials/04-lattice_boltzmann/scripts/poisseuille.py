## A script to simulate planar Poisseuille flow in Espresso
from espressomd import lb, shapes, lbboundaries
import numpy as np



# System setup
system = espressomd.System(box_l=[1.0, 1.0, 1.0])
system.seed  = system.cell_system.get_state()['n_nodes'] * [1234]
system.time_step = 0.01
system.cell_system.skin = 0.2

box_l = 16
system.box_l = [box_l] * 3

lbf = lb.LBFluidGPU(agrid=1, dens=1, visc=1, tau=0.01, ext_force=[0, 0.001, 0])
system.actors.add(lbf)
system.thermostat.set_lb(kT=0)

# Setup boundaries
walls = [lbboundaries.LBBoundary() for k in range(2)]
walls[0].set_params(shape=shapes.Wall(normal=[1,0,0], dist = 1.5))
walls[1].set_params(shape=shapes.Wall(normal=[-1,0,0], dist = -14.5))

for wall in walls:
    system.lbboundaries.add(wall)

## Perform enough iterations until the flow profile
## is static (5000 LB updates):
system.integrator.run(5000)

## Part of the solution
node_v_list = []
for i in range(box_l):
    node_v_list.append(lbf[i, 0, 0].velocity[1])

with open("lb_fluid_velocity.dat", "w") as f:
    for line in node_v_list:
        f.write(str(line)+"\n")
