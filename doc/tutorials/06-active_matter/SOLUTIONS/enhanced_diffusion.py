################################################################################
#                                                                              #
# Copyright (C) 2010-2018 The ESPResSo project                                 #
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
#                  Active Matter: Enhanced Diffusion Tutorial                  #
#                                                                              #
##########################################################################

from __future__ import print_function

import numpy as np
import os
import sys
import time

import espressomd
from espressomd import assert_features
from espressomd.observables import ParticlePositions, ParticleVelocities, ParticleAngularVelocities
from espressomd.accumulators import Correlator

required_features = ["ENGINE", "ROTATION"]
assert_features(required_features)

# create an output folder
outdir = "./RESULTS_ENHANCED_DIFFUSION/"
try:
    os.makedirs(outdir)
except:
    print("INFO: Directory \"{}\" exists".format(outdir))

##########################################################################

# Read in the active velocity from the command prompt
if len(sys.argv) != 2:
    print("Usage:", sys.argv[0], "<vel> (0 <= vel < 10.0)")
    exit()

vel = float(sys.argv[1])

# Set the basic simulation parameters
sampsteps = 5000
samplength = 1000
tstep = 0.01

system = espressomd.System(box_l=[10.0, 10.0, 10.0])
system.cell_system.skin = 0.3
system.time_step = tstep

################################################################################
#
# To obtain accurate statistics, you will need to run the simulation
# several times, which is accomplished by this loop. Do not increase
# this number too much, as it will slow down the simulation.
#
##########################################################################

for run in range(5):
    # Set up a random seed (a new one for each run)
    system.seed = np.random.randint(0, 2**31 - 1)

    # Use the Langevin thermostat (no hydrodynamics)
    system.thermostat.set_langevin(kT=1.0, gamma=1.0, seed=42)

    # Place a single active particle (that can rotate freely! rotation=[1,1,1])
    system.part.add(pos=[5.0, 5.0, 5.0], swimming={
                    'v_swim': vel}, rotation=[1, 1, 1])

    # Initialize the mean squared displacement (MSD) correlator
    tmax = tstep * sampsteps

    pos_id = ParticlePositions(ids=[0])
    msd = Correlator(obs1=pos_id,
                     corr_operation="square_distance_componentwise",
                     delta_N=1,
                     tau_max=tmax,
                     tau_lin=16)
    system.auto_update_accumulators.add(msd)

    # Initialize the velocity auto-correlation function (VACF) correlator
    vel_id = ParticleVelocities(ids=[0])
    vacf = Correlator(obs1=vel_id,
                      corr_operation="scalar_product",
                      delta_N=1,
                      tau_max=tmax,
                      tau_lin=16)
    system.auto_update_accumulators.add(vacf)

    # Initialize the angular velocity auto-correlation function (AVACF)
    # correlator
    ang_id = ParticleAngularVelocities(ids=[0])
    avacf = Correlator(obs1=ang_id,
                       corr_operation="scalar_product",
                       delta_N=1,
                       tau_max=tmax,
                       tau_lin=16)
    system.auto_update_accumulators.add(avacf)

    # Integrate 5,000,000 steps. This can be done in one go as well.
    for i in range(sampsteps):
        if (i + 1) % 100 == 0:
            print('\rrun %i: %.0f%%' % (run + 1, (i + 1) * 100. / sampsteps),
                  end='')
            sys.stdout.flush()
        system.integrator.run(samplength)
    print()

    # Finalize the accumulators and write to disk
    system.auto_update_accumulators.remove(msd)
    msd.finalize()
    np.savetxt("{}/msd_{}_{}.dat".format(outdir, vel, run), msd.result())

    system.auto_update_accumulators.remove(vacf)
    vacf.finalize()
    np.savetxt("{}/vacf_{}_{}.dat".format(outdir, vel, run), vacf.result())

    system.auto_update_accumulators.remove(avacf)
    avacf.finalize()
    np.savetxt("{}/avacf_{}_{}.dat".format(outdir, vel, run), avacf.result())
