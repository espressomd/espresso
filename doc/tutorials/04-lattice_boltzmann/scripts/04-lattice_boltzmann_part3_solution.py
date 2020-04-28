# Copyright (C) 2019 The ESPResSo project
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

import logging
import numpy as np
logging.basicConfig(level=logging.INFO)

import espressomd
espressomd.assert_features(['LENNARD_JONES'])
import espressomd.accumulators
import espressomd.observables
import espressomd.polymer
import espressomd.minimize_energy

# Setup constant
TIME_STEP = 0.01
LOOPS = 5000
STEPS = 100

# System setup
system = espressomd.System(box_l=[32.0, 32.0, 32.0])
system.cell_system.skin = 0.4

# Lennard-Jones interaction
system.non_bonded_inter[0, 0].lennard_jones.set_params(
    epsilon=1.0, sigma=1.0,
    shift="auto", cutoff=2.0**(1.0 / 6.0))

# Fene interaction
fene = espressomd.interactions.FeneBond(k=7, d_r_max=2)
system.bonded_inter.add(fene)

N_MONOMERS = [10, 20, 40]

rh_results = np.zeros(len(N_MONOMERS))
diffusion_constant_results = np.zeros(len(N_MONOMERS))
for index, N in enumerate(N_MONOMERS):
    logging.info("Polymer size: {}".format(N))
    system.part.clear()
    system.thermostat.turn_off()
    system.actors.clear()
    system.auto_update_accumulators.clear()
    # Setup polymer of part_id 0 with fene bond
    positions = espressomd.polymer.linear_polymer_positions(n_polymers=1,
                                                            beads_per_chain=N,
                                                            bond_length=1, seed=5642,
                                                            min_distance=0.9)
    for i, pos in enumerate(positions[0]):
        pid = len(system.part)
        system.part.add(id=pid, pos=pos)
        if i > 0:
            system.part[pid].add_bond((fene, pid - 1))

    logging.info("Warming up the polymer chain.")
    system.time_step = 0.002
    espressomd.minimize_energy.steepest_descent(
        system,
        f_max=1.0,
        gamma=10,
        max_steps=2000,
        max_displacement=0.01)
    logging.info("Warmup finished.")

    logging.info("Equilibration.")
    system.time_step = TIME_STEP
    system.thermostat.set_langevin(kT=1.0, gamma=10, seed=42)
    system.integrator.run(50000)
    logging.info("Equilibration finished.")

    system.thermostat.turn_off()

    lbf = espressomd.lb.LBFluidGPU(
        kT=1,
        seed=123,
        agrid=1,
        dens=1,
        visc=5,
        tau=TIME_STEP)
    system.actors.add(lbf)
    system.thermostat.set_lb(LB_fluid=lbf, seed=142, gamma=5)

    logging.info("Warming up the system with LB fluid.")
    system.integrator.run(1000)
    logging.info("Warming up the system with LB fluid finished.")

    # configure correlator
    com_pos = espressomd.observables.ComPosition(ids=range(N))
    correlator = espressomd.accumulators.Correlator(
        obs1=com_pos, tau_lin=16, tau_max=LOOPS * STEPS, delta_N=1,
        corr_operation="square_distance_componentwise", compress1="discard1")
    system.auto_update_accumulators.add(correlator)

    logging.info("Sampling started.")
    rhs = np.zeros(LOOPS)
    for i in range(LOOPS):
        system.integrator.run(STEPS)
        rhs[i] = system.analysis.calc_rh(
            chain_start=0,
            number_of_chains=1,
            chain_length=N)[0]

    logging.info("Sampling finished.")

    correlator.finalize()
    corrdata = correlator.result()

    rh_results[index] = np.average(rhs)
    tau = corrdata[1:, 0]
    msd = corrdata[1:, 2] + corrdata[1:, 3] + corrdata[1:, 4]
    np.save('msd_{}'.format(N), np.c_[tau, msd])

np.save('rh.npy', rh_results)
