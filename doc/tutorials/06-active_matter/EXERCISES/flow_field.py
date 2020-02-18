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
#            Active Matter: Swimmer Flow Field Tutorial
#
##########################################################################

import os

import espressomd
from espressomd import assert_features, lb


assert_features(["ENGINE", "CUDA", "MASS", "ROTATION", "ROTATIONAL_INERTIA"])


## Exercise 1 ##
# Create a routine to read in the hydrodynamic type
# (pusher/puller) and position at which the particle
# is initiated, set the variables 'type' and 'pos' to
# these values, respectively.

...

mode = ...
pos = ...

##########################################################################

## Exercise 2 ##
# Create an output directory that is labeled according
# to the value of the type and position, use the parameter
# 'outdir' to store this path

outdir = ...
os.makedirs(outdir, exist_ok=True)

# System parameters

LENGTH = 25.0
PROD_STEPS = 1000
PROD_LENGTH = 50
TIME_STEP = 0.01

system = espressomd.System(box_l=[LENGTH, LENGTH, LENGTH])
system.cell_system.skin = 0.3
system.time_step = TIME_STEP
system.min_global_cut = 1.0

##########################################################################

# Set the position of the particle

## Exercise 3 ##
# Determine the initial position of the particle, which
# should be in the center of the box, and shifted by
# the value of 'pos' in the direction of the z-axis

x0 = ...
y0 = ...
z0 = ...

# Sphere size, mass, and moment of inertia, dipole force

sph_size = 0.5
sph_mass = 4.8
Ixyz = 4.8
force = 0.1

## Exercise 4 ##
# Why is the sphere size set to 0.5 (this value is
# an approximation for the real value)? What happens when you
# change the mass and rotational inertia? Why is the value of
# the force chosen to be low.

# Setup the particle

system.part.add(
    pos=[x0, y0, z0], type=0, mass=sph_mass, rinertia=[Ixyz, Ixyz, Ixyz],
    swimming={'f_swim': force, 'mode': mode, 'dipole_length': sph_size + 0.5})

## Exercise 5 ##
# Why is the dipole_length chosen in this way?
# What happens if you make the length go to zero?
# Why does this happen?

##########################################################################

# Setup the fluid (quiescent)

lbf = lb.LBFluidGPU(agrid=1.0, dens=1.0, visc=1.0,
                    tau=TIME_STEP)
## Exercise 6 ##
# Can the particle rotate in the flow field?
system.actors.add(lbf)
system.thermostat.set_lb(LB_fluid=lbf, gamma=20.0, seed=42)

##########################################################################

# Output the coordinates

with open("{}/trajectory.dat".format(outdir), 'w') as outfile:
    print("####################################################", file=outfile)
    print("#        time        position       velocity       #", file=outfile)
    print("####################################################", file=outfile)

    # Production run

    for k in range(PROD_STEPS):
        # Output quantities
        print("{time} {pos[0]} {pos[1]} {pos[2]} {vel[0]} {vel[1]} {vel[2]}"
              .format(time=system.time, pos=system.part[0].pos, vel=system.part[0].v),
              file=outfile)

        # Output 50 simulations
        if k % (PROD_STEPS / 50) == 0:
            num = k / (PROD_STEPS / 50)
            lbf.print_vtk_velocity("{}/lb_velocity_{}.vtk".format(outdir, num))
            system.part.writevtk(
                "{}/position_{}.vtk".format(outdir, num), types=[0])

        system.integrator.run(PROD_LENGTH)

## Exercise 7 ##
# Use the snapshots and paraview to visualize the final state.
# By appropriately choosing the initial position, you can ensure
# that the swimmer is in the center of the box. Explain why
# the flow lines look the way they do.
