import espressomd

from espressomd import lb
from espressomd import lbboundaries
from espressomd import shapes
from espressomd.shapes import Rhomboid
from espressomd.shapes import Cylinder
from espressomd.interactions import OifLocalForces
from espressomd.interactions import OifGlobalForces
from espressomd.interactions import OifOutDirection

import numpy as np
import os, sys

from object_in_fluid import *

system = espressomd.System()
system.time_step = 0.1
system.box_l = [10.0, 10.0, 10.0]
system.cell_system.skin = 0.2

# creating the template for OIF object
cell_type = OifCellType(nodes_file="input/sphere393nodes.dat", triangles_file="input/sphere393triangles.dat", system = system, ks=1.0, kb=1.0, kal=1.0, kag=0.1, kv=0.1, check_orientation=False, resize=(3.0,3.0,3.0))

# creating the OIF object
cell0 = OifCell(cell_type=cell_type, part_type=0, origin=[5.0,5.0,5.0])
# cell0.output_vtk_pos_folded(file_name="cell0_0.vtk")


# fluid
lbf = espressomd.lb.LBFluid(agrid = 1, dens = 1.0, visc = 1.5, tau = 0.1, fric = 1.5)
system.actors.add(lbf)

diameter_init = cell0.diameter()
print("initial diameter = " + str(diameter_init))

# OIF object is being stretched by factor 1.5
maxCycle = 500
for p in system.part:
    newx = p.pos[0] - 5.0
    newx = newx * 1.5
    newx = newx + 5.0
    newy = p.pos[1] - 5.0
    newy = newy * 1.5
    newy = newy + 5.0
    newz = p.pos[2] - 5.0
    newz = newz * 1.5
    newz = newz + 5.0
    p.pos = [newx, newy, newz]

diameter_stretched = cell0.diameter()
print("stretched diameter = " + str(diameter_stretched))

# main integration loop
# OIF object is let to relax into relaxed shape of the sphere
for i in range(1,maxCycle):
    system.integrator.run(steps=10)
#    cell0.output_vtk_pos_folded(file_name="cell0_" + str(i) + ".vtk")


# final diameter is measured and compared to initial diameter
diameter_final = cell0.diameter()
print("final diameter = " + str(diameter_final))

# if the relaxed diameter is the same as initial diameter (up to 0.1 percent), test passed.
if ((diameter_final - diameter_init)/diameter_init < 0.001 and (diameter_init - diameter_final)/diameter_init > -0.001):
    print("Test passed.")
else:
    print("Test failed.")

