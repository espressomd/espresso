import espressomd

from espressomd import lb
from espressomd import lbboundaries
from espressomd import shapes
from espressomd.shapes import Rhomboid
from espressomd.shapes import Cylinder
from espressomd.interactions import OifLocalForces
from espressomd.interactions import OifGlobalForces
from espressomd.interactions import OifOutDirection
# from espressomd.interactions import SoftSphereInteraction
# from espressomd.interactions import MembraneCollisionInteraction

import numpy as np
import os, sys

from oif_classes import OifCellType
from oif_classes import OifCell
from oif_classes import Edge
from oif_classes import Triangle
from oif_classes import PartPoint
from oif_classes import FixedPoint

from oif_utils import AngleBtwTriangles, Distance, AreaTriangle, GetNTriangle, KS, OutputVtkRhomboid, OutputVtkCylinder

simNo = sys.argv[1]
os.mkdir("output/sim" + str(simNo))

system = espressomd.System()
system.time_step = 0.1
system.box_l = [22.0, 14.0, 15.0]
system.cell_system.skin = 0.2

# creating the template for RBCs
cell_type = OifCellType(nodesfile="input/rbc374nodes.dat", trianglesfile="input/rbc374triangles.dat", system = system, ks=0.04, kb=0.016, kal=0.02, kag=0.9, kv=1.0, checkOrientation=False, stretch=[2.0,2.0,2.0])

# creating the RBCs
cell0 = OifCell(cellType=cell_type, partType=0, origin=[5.0,5.0,3.0])
cell1 = OifCell(cellType=cell_type, partType=1, origin=[5.0,5.0,7.0])

# fluid
lbf = espressomd.lb.LBFluid(agrid = 1, dens = 1.0, visc = 1.5, tau = 0.1, fric = 1.5 , ext_force = [0.002, 0.0, 0.0])
system.actors.add(lbf)

# creating boundaries and obstacles in the channel
# OutputVtk writes a file
# lbboundaries created boundaries for fluid
# contraints created aboundaries for the cells

boundaries = []

boxX = 22.0
boxY = 14.0
boxZ = 15.0
time_step = 0.1

# bottom of the channel
boundaries.append(shapes.Rhomboid(corner=[0.0,0.0,0.0], a=[boxX,0.0,0.0], b=[0.0,boxY,0.0], c=[0.0,0.0,1.0], direction = 1))
OutputVtkRhomboid(corner=[0.0,0.0,0.0], a=[boxX,0.0,0.0], b=[0.0,boxY,0.0], c=[0.0,0.0,1.0], outFile="output/sim"+str(simNo)+"/wallBottom.vtk")

# top of the channel
boundaries.append(shapes.Rhomboid(corner=[0.0,0.0,boxZ-1], a=[boxX,0.0,0.0], b=[0.0,boxY,0.0], c=[0.0,0.0,1.0], direction = 1))
OutputVtkRhomboid(corner=[0.0,0.0,boxZ-1], a=[boxX,0.0,0.0], b=[0.0,boxY,0.0], c=[0.0,0.0,1.0], outFile="output/sim"+str(simNo)+"/wallTop.vtk")

# front wall of the channel
boundaries.append(shapes.Rhomboid(corner=[0.0,0.0,0.0], a=[boxX,0.0,0.0], b=[0.0,1.0,0.0], c=[0.0,0.0,boxZ], direction = 1))
OutputVtkRhomboid(corner=[0.0,0.0,0.0], a=[boxX,0.0,0.0], b=[0.0,1.0,0.0], c=[0.0,0.0,boxZ], outFile="output/sim"+str(simNo)+"/wallFront.vtk")

# back wall of the channel
boundaries.append(shapes.Rhomboid(corner=[0.0,boxY-1.0,0.0], a=[boxX,0.0,0.0], b=[0.0,1.0,0.0], c=[0.0,0.0,boxZ], direction = 1))
OutputVtkRhomboid(corner=[0.0,boxY-1.0,0.0], a=[boxX,0.0,0.0], b=[0.0,1.0,0.0], c=[0.0,0.0,boxZ], outFile="output/sim"+str(simNo)+"/wallBack.vtk")

# obstacle - cylinder A
boundaries.append(shapes.Cylinder(center=[11.0,2.0,7.0], axis=[0.0,0.0,1.0], length=7.0, radius=2.0, direction = 1))
OutputVtkCylinder(center=[11.0,2.0,7.0], axis=[0.0,0.0,1.0], length=7.0, radius=2.0, n=20, outFile="output/sim"+str(simNo)+"/cylinderA.vtk")

# obstacle - cylinder B
boundaries.append(shapes.Cylinder(center=[16.0,8.0,7.0], axis=[0.0,0.0,1.0], length=7.0, radius=2.0, direction = 1))
OutputVtkCylinder(center=[16.0,8.0,7.0], axis=[0.0,0.0,1.0], length=7.0, radius=2.0, n=20, outFile="output/sim"+str(simNo)+"/cylinderB.vtk")

# obstacle - cylinder C
boundaries.append(shapes.Cylinder(center=[11.0,12.0,7.0], axis=[0.0,0.0,1.0], length=7.0, radius=2.0, direction = 1))
OutputVtkCylinder(center=[11.0,12.0,7.0], axis=[0.0,0.0,1.0], length=7.0, radius=2.0, n=20, outFile="output/sim"+str(simNo)+"/cylinderC.vtk")

for boundary in boundaries:
    system.lbboundaries.add(lbboundaries.LBBoundary(shape = boundary))
    system.constraints.add(shape = boundary, particle_type = 10)

# cell-wall interactions
system.non_bonded_inter[0,10].soft_sphere.set_params(soft_a = 0.0001, soft_n = 1.2, soft_cut = 0.1, soft_offset = 0.0)
system.non_bonded_inter[1,10].soft_sphere.set_params(soft_a = 0.0001, soft_n = 1.2, soft_cut = 0.1, soft_offset = 0.0)

# cell-cell interactions
system.non_bonded_inter[0,1].membrane_collision.set_params(membrane_a = 0.0001, membrane_n = 1.2, membrane_cut = 0.1, membrane_offset = 0.0)

maxCycle = 50
# main integration loop
cell0.OutputVtkPosFolded(filename="output/sim" + str(simNo) + "/cell0_0.vtk")
cell1.OutputVtkPosFolded(filename="output/sim" + str(simNo) + "/cell1_0.vtk")
for i in range(1,maxCycle):
    system.integrator.run(steps=500)
    cell0.OutputVtkPosFolded(filename="output/sim" + str(simNo) + "/cell0_" + str(i) + ".vtk")
    cell1.OutputVtkPosFolded(filename="output/sim" + str(simNo) + "/cell1_" + str(i) + ".vtk")
    print "time: ", str(i*time_step)
print "Simulation completed."
