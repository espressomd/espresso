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
Simulate the motion of flexible red blood cells in a lattice-Boltzmann fluid
with solid obstacles. For more details, see :ref:`Object-in-fluid`.
"""

import espressomd

required_features = ["LB_BOUNDARIES", "EXTERNAL_FORCES", "SOFT_SPHERE",
                     "MASS"]
espressomd.assert_features(required_features)

from espressomd import lbboundaries
from espressomd import shapes

import os
import argparse
import warnings

import object_in_fluid as oif
from object_in_fluid.oif_utils import output_vtk_rhomboid, output_vtk_cylinder

parser = argparse.ArgumentParser()
parser.add_argument("sim", metavar="N", type=int, help="simulation identifier")
args = parser.parse_args()

output_path = "output/sim" + str(args.sim)
if not os.path.isdir(output_path):
    os.makedirs(output_path)
    print('Saving data to ' + output_path)
else:
    warnings.warn("Folder {} already exists, files will be overwritten"
                  .format(output_path))

boxX = 22.0
boxY = 14.0
boxZ = 15.0
time_step = 0.1

system = espressomd.System(box_l=(boxX, boxY, boxZ))
system.time_step = time_step
system.cell_system.skin = 0.2

# creating the template for RBCs
cell_type = oif.OifCellType(
    nodes_file="input/rbc374nodes.dat", triangles_file="input/rbc374triangles.dat",
    system=system, ks=0.04, kb=0.016, kal=0.02, kag=0.9, kv=1.0, check_orientation=False, resize=(2.0, 2.0, 2.0))

# creating the RBCs
cell0 = oif.OifCell(cell_type=cell_type,
                    particle_type=0, origin=[5.0, 5.0, 3.0])
cell1 = oif.OifCell(cell_type=cell_type,
                    particle_type=1, origin=[5.0, 5.0, 7.0])

# cell-wall interactions
system.non_bonded_inter[0, 10].soft_sphere.set_params(
    a=0.0001, n=1.2, cutoff=0.1, offset=0.0)
system.non_bonded_inter[1, 10].soft_sphere.set_params(
    a=0.0001, n=1.2, cutoff=0.1, offset=0.0)

# fluid
lbf = espressomd.lb.LBFluid(agrid=1, dens=1.0, visc=1.5, tau=0.1,
                            ext_force_density=[0.002, 0.0, 0.0])
system.actors.add(lbf)
system.thermostat.set_lb(LB_fluid=lbf, gamma=1.5)

# creating boundaries and obstacles in the channel
# OutputVtk writes a file
# lbboundaries created boundaries for fluid
# constraints created boundaries for the cells

boundaries = []

# bottom of the channel
bottom_shape = shapes.Rhomboid(corner=[0.0, 0.0, 0.0], a=[boxX, 0.0, 0.0],
                               b=[0.0, boxY, 0.0], c=[0.0, 0.0, 1.0],
                               direction=1)
boundaries.append(bottom_shape)
output_vtk_rhomboid(
    bottom_shape, out_file=output_path + "/wallBottom.vtk")

# top of the channel
top_shape = shapes.Rhomboid(corner=[0.0, 0.0, boxZ - 1], a=[boxX, 0.0, 0.0],
                            b=[0.0, boxY, 0.0], c=[0.0, 0.0, 1.0], direction=1)
boundaries.append(top_shape)
output_vtk_rhomboid(
    top_shape, out_file=output_path + "/wallTop.vtk")

# front wall of the channel
front_shape = shapes.Rhomboid(corner=[0.0, 0.0, 0.0], a=[boxX, 0.0, 0.0],
                              b=[0.0, 1.0, 0.0], c=[0.0, 0.0, boxZ], direction=1)
boundaries.append(front_shape)
output_vtk_rhomboid(
    front_shape, out_file=output_path + "/wallFront.vtk")

# back wall of the channel
back_shape = shapes.Rhomboid(corner=[0.0, boxY - 1.0, 0.0], a=[boxX, 0.0, 0.0],
                             b=[0.0, 1.0, 0.0], c=[0.0, 0.0, boxZ], direction=1)
boundaries.append(back_shape)
output_vtk_rhomboid(
    back_shape, out_file=output_path + "/wallBack.vtk")

# obstacle - cylinder A
cylA_shape = shapes.Cylinder(center=[11.0, 2.0, boxZ / 2.], axis=[0.0, 0.0, 1.0],
                             length=boxZ, radius=2.0, direction=1)
boundaries.append(cylA_shape)
output_vtk_cylinder(
    cylA_shape, n=20, out_file=output_path + "/cylinderA.vtk")

# obstacle - cylinder B
cylB_shape = shapes.Cylinder(center=[16.0, 8.0, boxZ / 2.], axis=[0.0, 0.0, 1.0],
                             length=boxZ, radius=2.0, direction=1)
boundaries.append(cylB_shape)
output_vtk_cylinder(
    cylB_shape, n=20, out_file=output_path + "/cylinderB.vtk")

# obstacle - cylinder C
cylC_shape = shapes.Cylinder(center=[11.0, 12.0, boxZ / 2.], axis=[0.0, 0.0, 1.0],
                             length=boxZ, radius=2.0, direction=1)
boundaries.append(cylC_shape)
output_vtk_cylinder(
    cylC_shape, n=20, out_file=output_path + "/cylinderC.vtk")

for boundary in boundaries:
    system.lbboundaries.add(lbboundaries.LBBoundary(shape=boundary))
    system.constraints.add(shape=boundary, particle_type=10)


maxCycle = 50
# main integration loop
cell0.output_vtk_pos_folded(
    file_name=output_path + "/cell0_0.vtk")
cell1.output_vtk_pos_folded(
    file_name=output_path + "/cell1_0.vtk")
for i in range(1, maxCycle):
    system.integrator.run(steps=500)
    cell0.output_vtk_pos_folded(
        file_name=output_path + "/cell0_" + str(i) + ".vtk")
    cell1.output_vtk_pos_folded(
        file_name=output_path + "/cell1_" + str(i) + ".vtk")
    print("time: {:.1f}".format(i * time_step))
print("Simulation completed.")
