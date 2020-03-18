#
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
#
##########################################################################
#
#            Active Matter: Rectification System Setup
#
##########################################################################

from math import cos, pi, sin
import numpy as np
import os

import espressomd
espressomd.assert_features(["CUDA", "LB_BOUNDARIES_GPU"])
from espressomd import lb
from espressomd.lbboundaries import LBBoundary
import espressomd.shapes
from espressomd.shapes import Cylinder


# Setup constants

outdir = "./RESULTS_RECTIFICATION"
os.makedirs(outdir, exist_ok=True)

# Setup the box (we pad the geometry to make sure
# the LB boundaries are away from the edges of the box)

LENGTH = 100
DIAMETER = 20
PADDING = 2
TIME_STEP = 0.01

# Setup the MD parameters

BOX_L = np.array(
    [LENGTH + 2 * PADDING,
     DIAMETER + 2 * PADDING,
     DIAMETER + 2 * PADDING])
system = espressomd.System(box_l=BOX_L)
system.cell_system.skin = 0.1
system.time_step = TIME_STEP
system.min_global_cut = 0.5

# Setup LB fluid

lbf = lb.LBFluidGPU(agrid=1.0, dens=1.0, visc=1.0, tau=TIME_STEP)
system.actors.add(lbf)

##########################################################################
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
        center=BOX_L / 2.,
        axis=[1, 0, 0], radius=DIAMETER / 2.0, length=LENGTH, direction=-1))
system.lbboundaries.add(cylinder)

# Setup walls

## Exercise 1 ##
# Set up two walls to cap the cylinder respecting the padding
# between them and the edge of the box, with a normal along the x-axis
wall = ...

# Setup cone

IRAD = 4.0
ANGLE = pi / 4.0
ORAD = (DIAMETER - IRAD) / sin(ANGLE)
SHIFT = 0.25 * ORAD * cos(ANGLE)

hollow_cone = LBBoundary(shape=espressomd.shapes.HollowConicalFrustum(center=[BOX_L[0] / 2.0 - 1.3 * SHIFT,
                                                                              BOX_L[1] /
                                                                              2.0,
                                                                              BOX_L[2] / 2.0],
                                                                      axis=[-1,
                                                                            0, 0],
                                                                      r1=ORAD,
                                                                      r2=IRAD,
                                                                      thickness=2.0,
                                                                      length=18,
                                                                      direction=1))
system.lbboundaries.add(hollow_cone)

##########################################################################

# Output the geometry

lbf.print_vtk_boundary("{}/boundary.vtk".format(outdir))

##########################################################################
## Exercise 2 ##
# Visualize this geometry using paraview

##########################################################################
