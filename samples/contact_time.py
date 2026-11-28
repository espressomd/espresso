#
# Copyright (C) 2025-2026 The ESPResSo project
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
Calculate the contact time between particle pairs in a Lennard-Jones fluid.
"""

import espressomd
import espressomd.accumulators
import espressomd.observables
import numpy as np

espressomd.assert_features(["LENNARD_JONES"])

# Simulation parameters
N_steps = 5000
seed = 42
contact_threshold = 1.
system = espressomd.System(box_l=[14., 14., 14.])
system.time_step = 0.01
system.cell_system.skin = 1.
rng = np.random.default_rng(seed)
system.part.add(pos=rng.random((10, 3)) * np.copy(system.box_l))
ids = np.copy(system.part.all().id)

system.non_bonded_inter[0, 0].lennard_jones.set_params(
    epsilon=1., sigma=1., shift="auto", cutoff=3.)

# relax the system
system.integrator.set_steepest_descent(f_max=0.05, gamma=1.,
                                       max_displacement=0.01)
system.integrator.run(1000)

# setup LD integration
system.thermostat.set_langevin(kT=1., gamma=1., seed=seed)
system.integrator.set_vv() 

# Setup observables to track pairwise distances and particle positions
pairwise_dist_obs = espressomd.observables.PairwiseDistances(ids=ids,
                                                             target_ids=ids)
particle_pos_obs = espressomd.observables.ParticlePositions(ids=ids)

# Setup the accumulators to track the contact times and the time series
contact_time_accumulator = espressomd.accumulators.ContactTimes(
    obs=pairwise_dist_obs, delta_N=1, contact_threshold=contact_threshold)

system.auto_update_accumulators.add(contact_time_accumulator)

# Run the simulation
system.integrator.run(N_steps)

# Get the contact time calculated on the fly
contact_times = contact_time_accumulator.contact_times()

print(f"The contact times are {contact_times}")
