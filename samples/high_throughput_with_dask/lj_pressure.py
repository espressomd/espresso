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
Obtain the average pressure of a Lennard-Jones liquid.
For use with ``dask_espresso``.
Adapted from :file:`samples/lj_liquid.py`.
"""

import espressomd
import numpy as np

import dask_espresso as de


# Note: Avoid print() in this script, as the standard output is used
# for data transfer to Dask. Use the logging module for status messages.
import logging
logger = logging.getLogger(__name__)


# Get parameters from Dask via the standard input stream
params = de.get_data_from_stdin()

logger.info("Parameters:", params)
n_particles = int(params["n_particles"])
volume_fraction = float(params["volume_fraction"])
n_steps = int(params["n_steps"])

required_features = ["LENNARD_JONES"]
espressomd.assert_features(required_features)

rng = np.random.default_rng()

# System
#############################################################

# Interaction parameters (Lennard-Jones)
#############################################################

lj_eps = 1.0  # LJ epsilon
lj_sig = 1.0  # particle diameter
lj_cut = lj_sig * 2**(1. / 6.)  # cutoff distance

# System parameters
#############################################################
# volume of N spheres with radius r: N * (4/3*pi*r^3)
box_l = (n_particles * 4. / 3. * np.pi * (lj_sig / 2.)**3
         / volume_fraction)**(1. / 3.)

# System
#############################################################
system = espressomd.System(box_l=3 * [box_l])

# Integration parameters
#############################################################
system.time_step = 0.01
system.cell_system.skin = 0.4

# Interaction setup
#############################################################
system.non_bonded_inter[0, 0].lennard_jones.set_params(
    epsilon=lj_eps, sigma=lj_sig, cutoff=lj_cut, shift="auto")

# Particle setup
#############################################################

system.part.add(pos=rng.random((n_particles, 3)) * system.box_l)

#  Warmup Integration
#############################################################

# warmup
logger.info("Warmup")
system.integrator.set_steepest_descent(
    f_max=0, max_displacement=0.01, gamma=1E-3)
system.integrator.run(1)
while np.any(np.abs(system.part.all().f) * system.time_step > .1): 
    system.integrator.run(10)

system.integrator.set_vv()
system.thermostat.set_langevin(kT=1.0, gamma=1.0, seed=42)

system.integrator.run(1000)
min_skin = 0.2
max_skin = 1.0
# tuning and equilibration
logger.info("Tune skin: {:.3f}".format(system.cell_system.tune_skin(
    min_skin=min_skin, max_skin=max_skin, tol=0.05, int_steps=100)))
system.integrator.run(1000)

logger.info("Measuring")

pressures = np.zeros(n_steps)
for i in range(n_steps):
    system.integrator.run(10)
    pressures[i] = system.analysis.pressure()["total"]

# Put the simulation results into a dictionary
result = {"pressure": np.mean(pressures),
          "pressure_std_dev": np.std(pressures)}

# Output the results for Dask via the standard output stream
print(de.encode_transport_data(result))
