## A script to simulate planar Poisseuille flow in Espresso
# a spherical RBC-like partice is added and advected
# with and without volume conservation
from espressomd import System, lb, shapes, lbboundaries
import numpy as np
from espressomd.virtual_sites import VirtualSitesInertialessTracers

# System setup
system = System()
system.time_step = 1/6.
system.cell_system.skin = 0.1
system.virtual_sites =VirtualSitesInertialessTracers()
print "Parallelization: " + str(system.cell_system.node_grid)
#system.cell_system.node_grid = [2, 2, 2]
#print system.cell_system.node_grid

boxZ = 20
system.box_l = [20, 20, boxZ]

force = 0.001
lbf = lb.LBFluid(agrid=1, dens=1, visc=1, tau= system.time_step, ext_force=[force, 0, 0], fric = 1)
system.actors.add(lbf)

system.thermostat.set_lb(kT=0)

# Setup boundaries
walls = [lbboundaries.LBBoundary() for k in range(2)]
walls[0].set_params(shape=shapes.Wall(normal=[0,0,1], dist = 0.5))
walls[1].set_params(shape=shapes.Wall(normal=[0,0,-1], dist = -boxZ + 0.5))

for wall in walls:
    system.lbboundaries.add(wall)
    
    
   
from addSoft import AddSoft
k1 = 0.1
k2 = 1
AddSoft(system, 10, 10, 10, k1, k2)

## case without bending and volCons
#outputDir = "outputPure"

## case with bending
from addBending import AddBending
kb = 1
AddBending(system, kb)
#outputDir = "outputBendPara"

## case with bending and volCons
from addVolCons import AddVolCons
kV = 10
AddVolCons(system, kV)
outputDir = "outputVolParaCUDA"

## make directory
from os import mkdir
mkdir(outputDir)
    
## Perform integration
from writeVTK import WriteVTK
WriteVTK(system, str(outputDir + "/cell_" + str(0) + ".vtk"))

stepSize = 1000
numSteps = 20

for i in range(0, numSteps):

    system.integrator.run(stepSize)
    WriteVTK(system, str(outputDir + "/cell_" + str(i+1) + ".vtk"))
    print "Done " + str(i+1) + " out of " + str(numSteps) + " steps."
