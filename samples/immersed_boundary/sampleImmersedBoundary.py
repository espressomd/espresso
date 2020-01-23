# Copyright (C) 2010-2019 The ESPResSo project
#
# This file is part of ESPResSo.
#
# ESPResSo is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ESPResSo is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""
Simulate the motion of a spherical red blood cell-like particle advected
in a planar Poiseuille flow, with or without volume conservation. For more
details, see :ref:`Immersed Boundary Method for soft elastic objects`.
"""
import espressomd

required_features = ["LB_BOUNDARIES", "VIRTUAL_SITES_INERTIALESS_TRACERS",
                     "EXPERIMENTAL_FEATURES"]
espressomd.assert_features(required_features)

from espressomd import lb, shapes, lbboundaries
from espressomd.virtual_sites import VirtualSitesInertialessTracers
import argparse

parser = argparse.ArgumentParser()
parser.add_argument(
    "--no-volcons", action="store_const", dest="volcons",
    const=False, help="Disable volume conservation", default=True)
parser.add_argument(
    "--no-bending", action="store_const", dest="bending",
    const=False, help="Disable bending", default=True)
args = parser.parse_args()
if args.volcons and not args.bending:
    print('Note: removing bending will also remove volume conservation')
    args.volcons = False

# System setup
boxZ = 20
system = espressomd.System(box_l=(20, 20, boxZ))
system.time_step = 1 / 6.
system.cell_system.skin = 0.1
system.virtual_sites = VirtualSitesInertialessTracers()
print("Parallelization: " + str(system.cell_system.node_grid))

force = 0.001
from addSoft import AddSoft
k1 = 0.1
k2 = 1
AddSoft(system, 10, 10, 10, k1, k2)

# case without bending and volCons
outputDir = "outputPure"

# case with bending
if args.bending:
    from addBending import AddBending
    kb = 1
    AddBending(system, kb)
    outputDir = "outputBendPara"

# case with bending and volCons
if args.volcons:
    from addVolCons import AddVolCons
    kV = 10
    AddVolCons(system, kV)
    outputDir = "outputVolParaCUDA"

# Add LB Fluid
lbf = lb.LBFluid(agrid=1, dens=1, visc=1, tau=system.time_step,
                 ext_force_density=[force, 0, 0])
system.actors.add(lbf)

system.thermostat.set_lb(LB_fluid=lbf, gamma=1.0, act_on_virtual=False)

# Setup boundaries
walls = [lbboundaries.LBBoundary() for k in range(2)]
walls[0].set_params(shape=shapes.Wall(normal=[0, 0, 1], dist=0.5))
walls[1].set_params(shape=shapes.Wall(normal=[0, 0, -1], dist=-boxZ + 0.5))

for wall in walls:
    system.lbboundaries.add(wall)

# make directory
import os
os.makedirs(outputDir)
print('Saving data to ' + outputDir)

# Perform integration
from writeVTK import WriteVTK
WriteVTK(system, str(outputDir + "/cell_" + str(0) + ".vtk"))

stepSize = 1000
numSteps = 20

for i in range(numSteps):
    system.integrator.run(stepSize)
    WriteVTK(system, str(outputDir + "/cell_" + str(i + 1) + ".vtk"))
    print("Done " + str(i + 1) + " out of " + str(numSteps) + " steps.")
