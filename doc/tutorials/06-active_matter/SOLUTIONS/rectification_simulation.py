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
#                     Active Matter: Rectification Tutorial                    #
#                                                                              #
##########################################################################

from __future__ import print_function

from math import cos, pi, sin
import numpy as np
import os
import sys

import espressomd
from espressomd import assert_features
from espressomd.shapes import Cylinder, Wall, HollowCone


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

# Read in the active velocity from the command prompt

if len(sys.argv) != 2:
    print("Usage:", sys.argv[0], "<vel> (0 <= vel < 10.0)")
    exit()

vel = float(sys.argv[1])

##########################################################################

# create an output folder

outdir = "./RESULTS_RECTIFICATION"
try:
    os.makedirs(outdir)
except:
    print("INFO: Directory \"{}\" exists".format(outdir))

# Setup the box (we pad the diameter to ensure that the LB boundaries
# and therefore the constraints, are away from the edge of the box)

length = 100
diameter = 20
prod_steps = 500
prod_length = 500
dt = 0.01

# Setup the MD parameters

system = espressomd.System(box_l=[length, diameter + 4, diameter + 4])
system.set_random_state_PRNG()
system.box_l = [length, diameter + 4, diameter + 4]
system.cell_system.skin = 0.1
system.time_step = dt
system.min_global_cut = 0.5
system.thermostat.set_langevin(kT=1.0, gamma=1.0, seed=42)

################################################################################
#
# Here we use exactly the same parameters for the geometry of the constraints
# that was used for the LB boundaries. This can be done, since the distance
# function used for the constraints is the same as the one used for the
# LB boundaries.
#
##########################################################################

cylinder = Cylinder(
    center=[length / 2.0, (diameter + 4) / 2.0, (diameter + 4) / 2.0],
    axis=[1, 0, 0], radius=diameter / 2.0, length=length, direction=-1)
system.constraints.add(shape=cylinder, particle_type=1)

# Setup walls

wall = Wall(dist=2, normal=[1, 0, 0])
system.constraints.add(shape=wall, particle_type=2)

wall = Wall(dist=-(length - 2), normal=[-1, 0, 0])
system.constraints.add(shape=wall, particle_type=3)

# Setup cone

irad = 4.0
angle = pi / 4.0
orad = (diameter - irad) / sin(angle)
shift = 0.25 * orad * cos(angle)

hollow_cone = HollowCone(
    center=[length / 2.0 + shift, (diameter + 4) / 2.0, (diameter + 4) / 2.0],
    axis=[-1, 0, 0],
    outer_radius=orad,
    inner_radius=irad,
    width=2.0,
    opening_angle=angle,
    direction=1)
system.constraints.add(shape=hollow_cone, particle_type=4)

################################################################################
#
# We set up a WCA (almost-hard) interaction between the particles and the
# confining geometry. We do not have particle-particle interactions, which
# are not necessary to observe rectification.
#
##########################################################################

sig = 0.5
cut = 1.12246 * sig
eps = 1.0
shift = 0.25

system.non_bonded_inter[0, 1].lennard_jones.set_params(
    epsilon=eps, sigma=sig, cutoff=cut, shift=shift)
system.non_bonded_inter[0, 2].lennard_jones.set_params(
    epsilon=eps, sigma=sig, cutoff=cut, shift=shift)
system.non_bonded_inter[0, 3].lennard_jones.set_params(
    epsilon=eps, sigma=sig, cutoff=cut, shift=shift)
system.non_bonded_inter[0, 4].lennard_jones.set_params(
    epsilon=eps, sigma=sig, cutoff=cut, shift=shift)

################################################################################
#
# Setup the particles. We put them all in two points one in each chamber
# and give them random directions. This speeds up the equilibration, since
# putting them all in a single chamber, would make it take a long time to
# observe the effect of rectification. Note that they need to be able to
# rotate freely, hence the command rotation=[1,1,1] is provided
#
##########################################################################

npart = 500
for cntr in range(npart):
    if cntr % 2 == 0:
        x = 0.25 * length
    else:
        x = 0.75 * length

    y = (diameter + 4) / 2.0
    z = (diameter + 4) / 2.0

    theta = float(2 * np.random.random() * np.pi)
    phi = float(2 * np.random.random() * np.pi)
    quats = a2quat(theta, phi)

    system.part.add(pos=[x, y, z], type=0, swimming={'v_swim': vel},
                    quat=quats, rotation=[1, 1, 1])

##########################################################################

# Equilibrate

system.integrator.run(25 * prod_length)

# Output the CMS coordinates

with open("{}/CMS_{}.dat".format(outdir, vel), "w") as outfile:
    print("####################################################", file=outfile)
    print("#        time    CMS x coord    average CMS        #", file=outfile)
    print("####################################################", file=outfile)

    # Production run

    dev_sum = 0.0
    dev_av = 0.0
    time_0 = system.time
    for i in range(prod_steps):
        # We output the coordinate of the center of mass in
        # the direction of the long axis, here we consider
        # the deviation from the center

        dev = system.galilei.system_CMS()[0] - 0.5 * length

        if i > 0:
            dev_sum = dev_sum + dev
            dev_av = dev_sum / i

        time = system.time - time_0

        print("{} {} {}".format(time, dev, dev_av), file=outfile)

        system.integrator.run(prod_length)

# Output the final configuration

system.part.writevtk("{}/points_{}.vtk".format(outdir, vel), types=[0])
