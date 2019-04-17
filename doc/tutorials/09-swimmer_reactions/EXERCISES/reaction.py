################################################################################
#                                                                              #
# Copyright (C) 2010-2018 The ESPResSo project        #
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
#               Catalytic Reactions: Enhanced Diffusion Tutorial               #
#                                                                              #
##########################################################################

from __future__ import print_function

import numpy as np
import os
import sys
import time
import espressomd
from espressomd import assert_features
assert_features(["ROTATION",
                 "ROTATIONAL_INERTIA",
                 "LANGEVIN_PER_PARTICLE",
                 "SWIMMER_REACTIONS",
                 "LENNARD_JONES"])
from espressomd.observables import ParticlePositions, ParticleBodyAngularVelocities
from espressomd.accumulators import Correlator
from espressomd.swimmer_reaction import Reaction

##########################################################################

# Read in the active velocity from the command prompt

if len(sys.argv) != 2:
    print("Usage:", sys.argv[0], "<passive/active = 0/1>")
    exit()

active = int(sys.argv[1])

if (active != 0) and (active != 1):
    print("Usage:", sys.argv[0], "<passive/active = 0/1>")
    exit()

# Set the parameters

box_l = 10
radius = 3.0
csmall = 0.1
rate = 1000.0

# Print input parameters

print("Box length: {}".format(box_l))
print("Colloid radius: {}".format(radius))
print("Particle concentration: {}".format(csmall))
print("Reaction rate: {}".format(rate))
print("Active or Passive: {}".format(active))

# Create output directory

if active == 0:
    outdir = "./passive-system"
else:
    outdir = "./active-system"

try:
    os.makedirs(outdir)
except:
    print("INFO: Directory \"{}\" exists".format(outdir))

##########################################################################

# Setup system parameters

equi_steps = 250
equi_length = 100

prod_steps = 2000
prod_length = 100

dt = 0.01

system = espressomd.System(box_l=[box_l, box_l, box_l])
system.cell_system.skin = 0.1
system.time_step = dt
system.min_global_cut = 1.1 * radius

# Set up the random seeds

system.set_random_state_PRNG()

##########################################################################

# Thermostat parameters

# Catalyzer is assumed to be larger, thus larger friction
frict_trans_colloid = 20.0
frict_rot_colloid = 20.0

# Particles are small and have smaller friction
frict_trans_part = 1.0
frict_rot_part = 1.0

# Temperature
temp = 1.0

##########################################################################

# Set up the swimmer

## Exercise 1 ##
# Determine the initial position of the particle, which
# should be in the center of the box.

x0pnt = ...
y0pnt = ...
z0pnt = ...

# Note that the swimmer needs to rotate freely
cent = len(system.part)
system.part.add(id=cent, pos=[x0pnt, y0pnt, z0pnt], type=0, temp=temp,
                gamma=frict_trans_colloid,
                gamma_rot=frict_rot_colloid,
                rotation=[1, 1, 1])

# Set up the particles

## Exercise 2 ##
# Above, we have set the concentration of the particles in the
# variable $csmall.  The concentration of both species of particles is
# equal.  Determine *how many* particles of one species there are.

# There are two species of equal concentration
nB = ...
nA = nB

print("Number of reactive A particles: {}".format(nB))
print("Number of reactive B particles: {}".format(nA))

for i in range(nA):
    x = box_l * np.random.random()
    y = box_l * np.random.random()
    z = box_l * np.random.random()

    # Prevent overlapping the colloid
    while (x - x0pnt)**2 + (y - y0pnt)**2 + (z - z0pnt)**2 < 1.15 * radius**2:
        x = box_l * np.random.random()
        y = box_l * np.random.random()
        z = box_l * np.random.random()

    # reactants and products do not need to rotate
    system.part.add(pos=[x, y, z], type=1, temp=temp,
                    gamma=frict_trans_part,
                    gamma_rot=frict_rot_part,
                    rotation=[0, 0, 0])

for i in range(nB):
    x = box_l * np.random.random()
    y = box_l * np.random.random()
    z = box_l * np.random.random()

    # Prevent overlapping the colloid
    while (x - x0pnt)**2 + (y - y0pnt)**2 + (z - z0pnt)**2 < 1.15 * radius**2:
        x = box_l * np.random.random()
        y = box_l * np.random.random()
        z = box_l * np.random.random()

    # reactants and products do not need to rotate
    system.part.add(pos=[x, y, z], type=2, temp=temp,
                    gamma=frict_trans_part,
                    gamma_rot=frict_rot_part,
                    rotation=[0, 0, 0])

print("box: {}, npart: {}".format(system.box_l, len(system.part)))

##########################################################################

# Set up the WCA potential

## Exercise 3 ##
# Why are there two different cutoff lengths for the LJ interaction
# catalyzer/product and catalyzer/reactant?


eps = 5.0
sig = 1.0
shift = 0.25
roff = radius - 0.5 * sig

# central and A particles
cut = 2**(1 / 6.) * sig
system.non_bonded_inter[0, 1].lennard_jones.set_params(
    epsilon=eps, sigma=sig, cutoff=cut, shift=shift, offset=roff)

# central and B particles (larger cutoff)
cut = 1.5 * sig
system.non_bonded_inter[0, 2].lennard_jones.set_params(
    epsilon=eps, sigma=sig, cutoff=cut, shift=shift, offset=roff)

##########################################################################

# Set up the reaction

cat_range = radius + 1.0 * sig
cat_rate = rate

## Exercise 4 ##
# We have read the activity parameter from the command line into
# $active, where 0 means off and 1 means on.  When $active = 0 we can
# simply go on, but when $active = 1 we have to set up the reaction.
# Check the $active parameter and setup a reaction for the catalyzer
# of type 0 with the reactants of type 1 and products of type 2.  The
# reaction range is stored in cat_range, the reaction rate in
# cat_rate.  Use the number-conserving scheme by setting swap on.

...
##########################################################################

# Perform warmup

cap = 20.0
warm_length = 100

## Exercise 5 ##
# Consult the User Guide for steepest descent integration to find out the
# difference to warmup with explicit force-capping.

system.integrator.set_steepest_descent(
    f_max=cap, gamma=.05, max_displacement=0.005)

# Note: system.part[:].f returns an n-by-3 array.
# We use the maximum force in any entry of this array as criterion
system.integrator.run(100)
while np.amax(system.part[:].f) >= cap:
    print(np.amax(system.part[:].f), np.argmax(
        system.part[:].f), cap, system.part[0].f)
    system.integrator.run(100)
# Restore the velocity Verlet integrator
system.integrator.set_vv()

##########################################################################

# Enable the thermostat

## Exercise 6 ##
# Why do we enable the thermostat only after warmup?

system.thermostat.set_langevin(kT=temp, gamma=frict_trans_colloid, seed=42)

##########################################################################

# Perform equilibration

# Integrate
for k in range(equi_steps):
    print("Equilibration: {} of {}".format(k, equi_steps))
    system.integrator.run(equi_length)

##########################################################################

for cnt in range(5):
    # Set up the MSD calculation

    tmax = prod_steps * prod_length * dt

    pos_id = ParticlePositions(ids=[cent])
    msd = Correlator(obs1=pos_id,
                     corr_operation="square_distance_componentwise",
                     delta_N=1,
                     tau_max=tmax,
                     tau_lin=16)
    system.auto_update_accumulators.add(msd)

    ## Exercise 7a ##
    # Construct the auto-correlators for the AVACF, using the example
    # of the MSD.
...

    # Perform production

    # Integrate
    for k in range(prod_steps):
        print("Production {} of 5: {} of {}".format(cnt + 1, k, prod_steps))
        system.integrator.run(prod_length)

    # Finalize the MSD and export
    system.auto_update_accumulators.remove(msd)
    msd.finalize()
    np.savetxt("{}/msd_{}.dat".format(outdir, cnt), msd.result())

    ## Exercise 7b ##
    # Finalize the angular velocity auto-correlation function (AVACF)
    # correlator and write the result to a file.
...
