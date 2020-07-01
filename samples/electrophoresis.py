#
# Copyright (C) 2013-2019 The ESPResSo project
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
Simulate electrophoresis of a linear polymer using the P3M electrostatics solver.
"""
import logging

import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize

import espressomd
import espressomd.observables
import espressomd.polymer
from espressomd import electrostatics, interactions

logging.basicConfig(level=logging.INFO)


required_features = ["P3M", "EXTERNAL_FORCES", "WCA"]
espressomd.assert_features(required_features)

N_SAMPLES = 1000
N_INT_STEPS = 100
E_FIELD = 1.0
N_MONOMERS = 20
N_IONS = 100
WARM_STEPS = 20
WARM_N_TIMES = 20
MIN_DIST = 0.9


system = espressomd.System(box_l=3 * [100.0])

system.time_step = 0.01
system.cell_system.skin = 0.4

# non-bonded interactions
###############################################################
# WCA between monomers
system.non_bonded_inter[0, 0].wca.set_params(epsilon=1, sigma=1)

# WCA counter-ions - polymer
system.non_bonded_inter[0, 1].wca.set_params(epsilon=1, sigma=1)

# WCA ions - polymer
system.non_bonded_inter[0, 2].wca.set_params(epsilon=1, sigma=1)

# WCA between ions
system.non_bonded_inter[1, 2].wca.set_params(epsilon=1, sigma=1)


# bonded interactions
################################################################
harmonic_bond = interactions.HarmonicBond(k=10, r_0=2)
angle_harmonic_bond = interactions.AngleHarmonic(bend=10, phi0=np.pi)
system.bonded_inter.add(harmonic_bond)
system.bonded_inter.add(angle_harmonic_bond)


# create monomer beads and bonds
##########################################################################
init_polymer_pos = espressomd.polymer.linear_polymer_positions(n_polymers=1, beads_per_chain=N_MONOMERS, bond_length=2.0,
                                                               seed=2, bond_angle=np.pi, min_distance=1.8, start_positions=np.array([system.box_l / 2.0]))

system.part.add(pos=init_polymer_pos[0], q=-np.ones(N_MONOMERS))

for i in range(1, N_MONOMERS):
    system.part[i].add_bond((harmonic_bond, i - 1))

for i in range(1, N_MONOMERS - 1):
    system.part[i].add_bond((angle_harmonic_bond, i - 1, i + 1))


# create counter-ions
###################################################################
system.part.add(pos=np.random.random((N_MONOMERS, 3)) * system.box_l,
                q=np.ones(N_MONOMERS),
                type=np.ones(N_MONOMERS, dtype=int))

# create excess ions
###############################################################
system.part.add(pos=np.random.random((N_IONS, 3)) * system.box_l,
                q=np.hstack((np.ones(N_IONS // 2), -np.ones(N_IONS // 2))),
                type=np.array(np.hstack((np.ones(N_IONS // 2), 2 * np.ones(N_IONS // 2))), dtype=int))

logging.info("particle types: {}\n".format(system.part[:].type))
logging.info("total charge: {}".format(np.sum(system.part[:].q)))

# warm-up integration
###############################################################
system.integrator.set_steepest_descent(f_max=0, gamma=1e-3,
                                       max_displacement=0.01)
i = 0
while system.analysis.min_dist() < MIN_DIST and i < WARM_N_TIMES:
    logging.debug(
        "total energy: {:+.2e}".format(system.analysis.energy()["total"]))
    system.integrator.run(WARM_STEPS)
    i += 1

logging.info(
    "total energy after warm-up: {:+.2e}\n".format(system.analysis.energy()["total"]))
system.integrator.set_vv()

system.thermostat.set_langevin(kT=1.0, gamma=1.0, seed=42)

# activate electrostatics
#############################################################
p3m = electrostatics.P3M(prefactor=1.0, accuracy=1e-2)
system.actors.add(p3m)

# apply external force (external electric field)
#############################################################
n_part = len(system.part)
system.part[:].ext_force = np.dstack(
    (system.part[:].q * np.ones(n_part) * E_FIELD, np.zeros(n_part), np.zeros(n_part)))[0]

# equilibration
#############################################################
system.integrator.run(500)

# observables for core analysis
#############################################################
obs_persistence_angles = espressomd.observables.CosPersistenceAngles(
    ids=system.part[:N_MONOMERS].id)
acc_persistence_angles = espressomd.accumulators.MeanVarianceCalculator(
    obs=obs_persistence_angles, delta_N=1)
system.auto_update_accumulators.add(acc_persistence_angles)

obs_bond_length = espressomd.observables.ParticleDistances(
    ids=system.part[:N_MONOMERS].id)
acc_bond_length = espressomd.accumulators.MeanVarianceCalculator(
    obs=obs_bond_length, delta_N=1)
system.auto_update_accumulators.add(acc_bond_length)

# data storage for python analysis
#############################################################
pos = np.full((N_SAMPLES, N_MONOMERS, 3), np.nan)

# sampling Loop
#############################################################
for i in range(N_SAMPLES):
    if i % 100 == 0:
        logging.info("\rsampling: {:4d}".format(i))
    system.integrator.run(N_INT_STEPS)
    pos[i] = system.part[:N_MONOMERS].pos

logging.info("\nsampling finished!\n")

# data analysis
############################################################

# calculate center of mass (COM) and its velocity
#############################################################
COM = pos.sum(axis=1) / N_MONOMERS
COM_v = (COM[1:] - COM[:-1]) / (N_INT_STEPS * system.time_step)

# calculate the electrophoretic mobility mu = v/E
##################################
mu = np.average(np.linalg.norm(COM_v, axis=1)) / E_FIELD
logging.info("electrophoretic mobility: {}".format(mu))

# calculate the persistence length...
#############################################################

# ...first python analysis
total_sampling_positions = []
total_cos_thetas = []
for positions in pos:
    bond_vectors = positions[1:, :] - positions[:-1, :]
    bond_lengths = np.linalg.norm(bond_vectors, axis=1)
    normed_bond_vectors = bond_vectors / bond_lengths[:, np.newaxis]
    # positions at which the angles between bonds are actually measured
    sampling_positions = np.insert(np.cumsum(bond_lengths)[:-1], 0, 0.0)
    cos_thetas = np.zeros_like(sampling_positions)
    for i in range(len(normed_bond_vectors)):
        cos_thetas[i] = np.dot(normed_bond_vectors[0], normed_bond_vectors[i])
    total_sampling_positions.append(sampling_positions)
    total_cos_thetas.append(cos_thetas)

sampling_positions = np.average(np.array(total_sampling_positions), axis=0)
cos_thetas = np.average(np.array(total_cos_thetas), axis=0)


def exponential(x, lp):
    return np.exp(-x / lp)


opt, _ = scipy.optimize.curve_fit(exponential, sampling_positions, cos_thetas)
persistence_length = opt[0]
logging.info("persistence length (python analysis): {}".format(
    persistence_length))

# ...second by using observables


def persistence_length_obs(
        acc_bond_length, acc_persistence_angles, exponential):
    bond_lengths_obs = np.array(acc_bond_length.get_mean())
    sampling_positions_obs = np.insert(
        np.cumsum(bond_lengths_obs)[:-1], 0, 0.0)
    cos_thetas_obs = np.array(acc_persistence_angles.get_mean())
    cos_thetas_obs = np.insert(cos_thetas_obs, 0, 1.0)

    opt_obs, _ = scipy.optimize.curve_fit(
        exponential, sampling_positions_obs, cos_thetas_obs)
    return sampling_positions_obs, cos_thetas_obs, opt_obs[0]


sampling_positions_obs, cos_thetas_obs, persistence_length_obs = persistence_length_obs(
    acc_bond_length, acc_persistence_angles, exponential)
logging.info("persistence length (observables): {}".format(
    persistence_length_obs))

# plot the results
#############################################################
fig, axs = plt.subplots(3)
axs[0].plot(COM[:, 0], label="COM pos in x-direction")
axs[0].plot(COM[:, 1], label="COM pos in y-direction")
axs[0].plot(COM[:, 2], label="COM pos in z-direction")
axs[0].legend()
axs[0].set_xlabel("time step")
axs[0].set_ylabel("r")

axs[1].plot(COM_v[:, 0], label="COM v in x-direction")
axs[1].plot(COM_v[:, 1], label="COM v in y-direction")
axs[1].plot(COM_v[:, 2], label="COM v in z-direction")
axs[1].legend()
axs[1].set_xlabel("time step")
axs[1].set_ylabel("v")

axs[2].plot(sampling_positions, cos_thetas, 'o',
            label="python analysis raw data")
axs[2].plot(sampling_positions_obs, cos_thetas_obs, 'o',
            label="observable analysis raw data")
axs[2].plot(sampling_positions, exponential(sampling_positions,
                                            opt[0]), label="exponential fit with python analysis")
axs[2].plot(sampling_positions_obs, exponential(sampling_positions_obs,
                                                persistence_length_obs), label="exponential fit with observable data")
axs[2].legend()
axs[2].set_xlabel("distance along polymer")
axs[2].set_ylabel(r"$\langle \cos(\theta) \rangle$")

plt.show()
