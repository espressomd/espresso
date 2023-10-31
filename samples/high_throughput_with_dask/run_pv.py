#
# Copyright (C) 2013-2023 The ESPResSo project
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

"""
Run a large amount of short ESPResSo simulations with the Dask parallel computing
library :cite:`rocklin15a`. Print the mean and standard error of the mean of the
scalar pressure of a Lennard-Jones gas simulated at different volume fractions.
"""

import sys
import dask.distributed
import logging
import numpy as np

import dask_espresso

logging.basicConfig(level=logging.WARN)

PYPRESSO = "/data/weeber/es/build/pypresso"  # adapt this path
SIM_SCRIPT = "lj_pressure.py"

# Simulation parameters
N_STEPS = int(2E4)
N_PARTICLES = 100
VOLUME_FRACTIONS = np.arange(0.1, 0.52, 0.01)


# Scheduler address
if len(sys.argv) != 2:
    raise Exception("Pass the scheduler address as command-line argument")
scheduler_address = sys.argv[1]

# Connect to scheduler
# Note: We pass a scheduler address here.
# dask.distributed.LocalCluster cannot be used, but clusters with
# remote workers such as HTCondorCluster likely can.
client = dask.distributed.Client(scheduler_address)


# List of futures for simulation results
futures = []

# Launch simulations asynchronously
for volume_fraction in VOLUME_FRACTIONS:
    sim_params = {"volume_fraction": volume_fraction,
                  "n_particles": N_PARTICLES,
                  "n_steps": N_STEPS}
    futures.append(client.compute(dask_espresso.dask_espresso_task(
        PYPRESSO, SIM_SCRIPT, **sim_params)))

# Show progress of calculation (optional)
dask.distributed.progress(futures)

# Gather the results of all futures, waiting for completion, if necessary
sim_results = client.gather(futures)

# Display results
for vol_frac, out in zip(VOLUME_FRACTIONS, sim_results):
    print(f"{vol_frac:3f} {out['pressure']:.3f} {out['pressure_std_dev']:.3f}")
