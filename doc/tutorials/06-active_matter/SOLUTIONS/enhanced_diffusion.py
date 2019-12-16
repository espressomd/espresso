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
#            Active Matter: Enhanced Diffusion Tutorial
#
##########################################################################

import numpy as np
import os
import argparse

import espressomd
from espressomd.observables import ParticlePositions, ParticleVelocities, ParticleAngularVelocities
from espressomd.accumulators import Correlator

espressomd.assert_features(
    ["ENGINE", "ROTATION", "MASS", "ROTATIONAL_INERTIA"])

# create an output folder

outdir = "./RESULTS_ENHANCED_DIFFUSION/"
os.makedirs(outdir, exist_ok=True)


##########################################################################

# Read in the active velocity from the command prompt

parser = argparse.ArgumentParser()
parser.add_argument('vel', type=float, help='Velocity of active particle.')
args = parser.parse_args()
vel = args.vel

# Set the basic simulation parameters

SAMP_STEPS = 5000
SAMP_LENGTH = 1000
T_STEP = 0.01

system = espressomd.System(box_l=[10.0, 10.0, 10.0])
system.cell_system.skin = 0.3
system.time_step = T_STEP

##########################################################################
#
# To obtain accurate statistics, you will need to run the simulation
# several times, which is accomplished by this loop. Do not increase
# this number too much, as it will slow down the simulation.
#
##########################################################################

for run in range(5):
    # Use the Langevin thermostat (no hydrodynamics)
    # Set up a random seed (a new one for each run)
    system.thermostat.set_langevin(kT=1.0, gamma=1.,
                                   seed=np.random.randint(0, 1000))

    # Place a single active particle that can rotate in all 3 dimensions.
    # Set mass and rotational inertia to separate the timescales for
    # translational and rotational diffusion.
    system.part.add(pos=[5.0, 5.0, 5.0], swimming={'v_swim': vel},
                    mass=0.1, rotation=3 * [True], rinertia=3 * [1.])

    # Initialize the mean squared displacement (MSD) correlator
    tmax = system.time_step * SAMP_STEPS

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
    for i in range(SAMP_STEPS):
        if (i + 1) % 100 == 0:
            print('\rrun %i: %.0f%%' % (run + 1, (i + 1) * 100. / SAMP_STEPS),
                  end='', flush=True)
        system.integrator.run(SAMP_LENGTH)
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
