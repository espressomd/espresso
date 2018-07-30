"""
This sample simulates planar Poisseuille flow in Espresso. A spherical RBC-like particle is added and advected with and without volume conservation.
"""
import espressomd

required_features = ["LB","LB_BOUNDARIES"]
espressomd.assert_features(required_features)

from espressomd import System, lb, shapes, lbboundaries
import numpy as np
from espressomd.virtual_sites import VirtualSitesInertialessTracers

# System setup
boxZ = 20
system = System(box_l=(20,20,boxZ))
system.time_step = 1/6.
system.cell_system.skin = 0.1
system.virtual_sites = VirtualSitesInertialessTracers()
print("Parallelization: " + str(system.cell_system.node_grid))

force = 0.001
lbf = lb.LBFluid(agrid=1, dens=1, visc=1, tau= system.time_step, ext_force_density=[force, 0, 0], fric = 1)
system.actors.add(lbf)

system.thermostat.set_lb(kT=0,act_on_virtual=False)

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
    print("Done " + str(i+1) + " out of " + str(numSteps) + " steps.")
