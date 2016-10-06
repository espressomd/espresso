## A script to simulate planar Poisseuille flow in Espresso
from espressomd import System, lb
import numpy as np

# System setup
system = System()
system.time_step = 0.1
system.cell_system.skin = 0.2
system.box_l = [16, 16, 16]

lbf = lb.LBFluid_GPU(agrid=1, dens=1, visc=1, tau=0.01, ext_force=[0, 0.001, 0])
system.actors.add(lbf)
system.thermostat.set_lb(kT=0)

# TODO
# ## Set LB boundaries
# lbboundary wall dist 1.5 normal 1. 0. 0.
# lbboundary wall dist -14.5 normal -1. 0. 0.

## Perform enough iterations until the flow profile
## is stable (1000 LB updates):
system.integrator.run(1e5)

## Part of the solution
node_v_list = []
for i in range(16):
    node_v_list.append(lbf.lbnode_get_node_velocity([i, 0, 0])[1])

np.savetxt("lb_fluid_velocity.dat", np.column_stack(range(16), node_v_list))
