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
"""
Set up a DPD fluid and calculate pressure as a function of the
varying density. The fluid is thermalized using a DPD thermostat.
"""

import espressomd

required_features = ["DPD"]
espressomd.assert_features(required_features)

import numpy as np

# Set up the box and time step
system = espressomd.System(box_l=3 * [10])
system.time_step = 0.01
system.cell_system.skin = 0.4

# DPD parameters
n_part = 200
kT = 1.
gamma = 1.5
r_cut = 1.
# Repulsive parameter
F_max = 1.

# Activate the thermostat
system.thermostat.set_dpd(kT=kT, seed=123)
np.random.seed(seed=42)

# Set up the DPD friction interaction
system.non_bonded_inter[0, 0].dpd.set_params(
    weight_function=0, gamma=gamma, r_cut=r_cut,
    trans_weight_function=0, trans_gamma=gamma, trans_r_cut=r_cut)

# Set up the repulsive interaction
system.non_bonded_inter[0, 0].hat.set_params(F_max=F_max, cutoff=r_cut)

# Add particles that are randomly distributed over the box
system.part.add(pos=system.box_l * np.random.random((n_part, 3)))

# As a usage example, we calculate the pressure at several
# particle densities.
sample_size = 100
int_steps = 1000
for V in range(100, 1000, 100):
    # Rescale the system to the new volume
    system.change_volume_and_rescale_particles(V**0.3333)

    # List of samples
    p_samples = []
    for i in range(sample_size):
        system.integrator.run(int_steps)
        p_samples.append(system.analysis.pressure()['total'])

    # Average pressure
    p_avg = np.mean(p_samples)
    # Standard deviation of pressure
    p_std = np.std(p_samples)

    print('rho {:.2f} p {:.2f} ({:.2f})'
          .format(float(n_part) / V, p_avg, p_std))
