################################################################################
#                                                                              #
# Copyright (C) 2010-2018 The ESPResSo project            #
#                                                                              #
# This file is part of ESPResSo.                                               #
#                                                                              #
# ESPResSo is free software: you can redistribute it and/or modify             #
# it under the terms of the GNU General Public License as published by         #
# the Free Software Foundation, either version 3 of the License, or            #
# (at your option) any later version.                                          #
#                                                                              #
# ESPResSo is distributed in the hope that it will be useful,                  #
# but WITHOUT ANY WARRANTY; without even the implied warranty of               #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                #
# GNU General Public License for more details.                                 #
#                                                                              #
# You should have received a copy of the GNU General Public License            #
# along with this program.  If not, see <http://www.gnu.org/licenses/>.        #
#                                                                              #
################################################################################
#                                                                              #
#                   Active Matter: Rectification System Setup                  #
#                                                                              #
##########################################################################

from __future__ import print_function

from math import cos, pi, sin
import numpy as np
import os
import sys

import espressomd
espressomd.assert_features(["CUDA", "LB_BOUNDARIES_GPU"])
from espressomd import assert_features, lb
from espressomd.lbboundaries import LBBoundary
from espressomd.shapes import Cylinder, Wall, HollowCone


# Setup constants

outdir = "./RESULTS_RECTIFICATION"
try:
    os.makedirs(outdir)
except:
    print("INFO: Directory \"{}\" exists".format(outdir))

# Setup the box (we pad the diameter to ensure that the LB boundaries
# and therefore the constraints, are away from the edge of the box)

length = 100
diameter = 20
dt = 0.01

# Setup the MD parameters

system = espressomd.System(box_l=[length, diameter + 4, diameter + 4])
system.cell_system.skin = 0.1
system.time_step = dt
system.min_global_cut = 0.5

# Setup LB parameters (these are irrelevant here) and fluid

agrid = 1
vskin = 0.1
frict = 20.0
visco = 1.0
densi = 1.0

lbf = lb.LBFluidGPU(agrid=agrid, dens=densi, visc=visco, tau=dt)
system.actors.add(lbf)
system.thermostat.set_lb(LB_fluid=lbf, gamma=frict, seed=42)

################################################################################
#
# Now we set up the three LB boundaries that form the rectifying geometry.
# The cylinder boundary/constraint is actually already capped, but we put
# in two planes for safety's sake. If you want to create a cylinder of
# 'infinite length' using the periodic boundaries, then the cylinder must
# extend over the boundary.
#
##########################################################################

# Setup cylinder

cylinder = LBBoundary(
    shape=Cylinder(
        center=[length / 2.0, (diameter + 4) / 2.0, (diameter + 4) / 2.0],
        axis=[1, 0, 0], radius=diameter / 2.0, length=length, direction=-1))
system.lbboundaries.add(cylinder)

# Setup walls

wall = LBBoundary(shape=Wall(dist=2, normal=[1, 0, 0]))
system.lbboundaries.add(wall)

wall = LBBoundary(shape=Wall(dist=-(length - 2), normal=[-1, 0, 0]))
system.lbboundaries.add(wall)

# Setup cone

irad = 4.0
angle = pi / 4.0
orad = (diameter - irad) / sin(angle)
shift = 0.25 * orad * cos(angle)

hollow_cone = LBBoundary(
    shape=HollowCone(
        center=[length / 2.0 + shift,
                (diameter + 4) / 2.0,
                (diameter + 4) / 2.0],
        axis=[-1, 0, 0],
        outer_radius=orad,
        inner_radius=irad,
        width=2.0,
        opening_angle=angle,
        direction=1))
system.lbboundaries.add(hollow_cone)

##########################################################################

# Output the geometry

lbf.print_vtk_boundary("{}/boundary.vtk".format(outdir))

##########################################################################
