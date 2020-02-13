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
#            Active Matter: Rectification Tutorial
#
##########################################################################

from math import cos, pi, sin
import numpy as np
import os
import argparse

import espressomd
from espressomd import assert_features
import espressomd.shapes
from espressomd.shapes import Cylinder, Wall


assert_features(["ENGINE", "LENNARD_JONES", "ROTATION", "MASS"])


# Quaternion procedure

def a2quat(phi, theta):

    q1w = cos(theta / 2.0)
    q1x = 0
    q1y = sin(theta / 2.0)
    q1z = 0

    q2w = cos(phi / 2.0)
    q2x = 0
    q2y = 0
    q2z = sin(phi / 2.0)

    q3w = (q1w * q2w - q1x * q2x - q1y * q2y - q1z * q2z)
    q3x = (q1w * q2x + q1x * q2w - q1y * q2z + q1z * q2y)
    q3y = (q1w * q2y + q1x * q2z + q1y * q2w - q1z * q2x)
    q3z = (q1w * q2z - q1x * q2y + q1y * q2x + q1z * q2w)

    return [q3w, q3x, q3y, q3z]


##########################################################################

parser = argparse.ArgumentParser()
parser.add_argument("vel", type=float, help="Velocity of active particles.")
args = parser.parse_args()

vel = args.vel
##########################################################################

# create an output folder
outdir = "./RESULTS_RECTIFICATION"
os.makedirs(outdir, exist_ok=True)

# Setup the box (we pad the diameter to ensure that the constraints
# are away from the edge of the box)

LENGTH = 100
DIAMETER = 20
PADDING = 2
PROD_STEPS = 500
PROD_LENGTH = 500
TIME_STEP = 0.005

# Setup the MD parameters

BOX_L = np.array(
    [LENGTH + 2 * PADDING,
     DIAMETER + 2 * PADDING,
     DIAMETER + 2 * PADDING])
system = espressomd.System(box_l=BOX_L)
system.cell_system.skin = 0.1
system.time_step = TIME_STEP
system.min_global_cut = 0.5
system.thermostat.set_langevin(kT=1.0, gamma=1.0, seed=42)

##########################################################################
#
# Here we use exactly the same parameters for the geometry of the constraints
# that was used for the LB boundaries. This can be done, since the distance
# function used for the constraints is the same as the one used for the
# LB boundaries.
#
##########################################################################

cylinder = Cylinder(
    center=0.5 * BOX_L,
    axis=[1, 0, 0], radius=DIAMETER / 2.0, length=LENGTH, direction=-1)
system.constraints.add(shape=cylinder, particle_type=1)

# Setup walls
wall = Wall(dist=PADDING, normal=[1, 0, 0])
system.constraints.add(shape=wall, particle_type=1)

wall = Wall(dist=-(LENGTH + PADDING), normal=[-1, 0, 0])
system.constraints.add(shape=wall, particle_type=1)

# Setup cone
IRAD = 4.0
ANGLE = pi / 4.0
ORAD = (DIAMETER - IRAD) / sin(ANGLE)
SHIFT = 0.25 * ORAD * cos(ANGLE)

conical_shape = espressomd.shapes.HollowConicalFrustum(center=[BOX_L[0] / 2.0 - 1.3 * SHIFT,
                                                               BOX_L[1] / 2.0,
                                                               BOX_L[2] / 2.0],
                                                       axis=[-1, 0, 0],
                                                       r1=ORAD,
                                                       r2=IRAD,
                                                       thickness=2.0,
                                                       length=18,
                                                       direction=1)
system.constraints.add(shape=conical_shape, particle_type=1)
##########################################################################
#
# We set up a WCA (almost-hard) interaction between the particles and the
# confining geometry. We do not have particle-particle interactions, which
# are not necessary to observe rectification.
#
##########################################################################

SIGMA = 0.5
CUTOFF = 1.12246 * SIGMA
EPSILON = 1
SHIFT = 0.25

system.non_bonded_inter[0, 1].lennard_jones.set_params(
    epsilon=EPSILON, sigma=SIGMA, cutoff=CUTOFF, shift=SHIFT)

##########################################################################
#
# Setup the particles. We put them all in two points one in each chamber
# and give them random directions. This speeds up the equilibration, since
# putting them all in a single chamber, would make it take a long time to
# observe the effect of rectification. Note that they need to be able to
# rotate freely, hence the command rotation=[1,1,1] is provided
#
##########################################################################

N_PART = 500
for cntr in range(N_PART):
    if cntr % 2 == 0:
        x = 0.25 * BOX_L[0]
    else:
        x = 0.75 * BOX_L[0]

    y = BOX_L[1] / 2.
    z = BOX_L[2] / 2.

    theta = float(np.random.random() * np.pi)
    phi = float(2 * np.random.random() * np.pi)
    quats = a2quat(theta, phi)

    system.part.add(pos=[x, y, z], type=0, swimming={'v_swim': vel},
                    quat=quats, rotation=[1, 1, 1])

##########################################################################
# Equilibrate
system.integrator.run(25 * PROD_LENGTH)

# Output the CMS coordinates
with open("{}/CMS_{}.dat".format(outdir, vel), "w") as outfile:
    print("####################################################", file=outfile)
    print("#        time    CMS x coord    average CMS        #", file=outfile)
    print("####################################################", file=outfile)

    # Production run
    dev_sum = 0.0
    dev_av = 0.0
    system.time = 0.
    for i in range(PROD_STEPS):
        if (i + 1) % 5 == 0:
            print('\rprogress: %.0f%%' % ((i + 1) * 100. / PROD_STEPS),
                  end='', flush=True)

        # We output the coordinate of the center of mass in
        # the direction of the long axis, here we consider
        # the deviation from the center (keep the padding in mind)

        dev = system.galilei.system_CMS()[0] - 0.5 * BOX_L[0]

        if i > 0:
            dev_sum = dev_sum + dev
            dev_av = dev_sum / i

        print("{} {} {}".format(system.time, dev, dev_av), file=outfile)

        system.integrator.run(PROD_LENGTH)

print()

# Output the final configuration

system.part.writevtk("{}/points_{}.vtk".format(outdir, vel), types=[0])
