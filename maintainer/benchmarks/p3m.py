#
# Copyright (C) 2013-2018 The ESPResSo project
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
from __future__ import print_function
import os
import sys
import numpy as np
import espressomd
from time import time, sleep
from espressomd import thermostat
from espressomd import electrostatics

required_features = ["ELECTROSTATICS", "LENNARD_JONES", "MASS"]
espressomd.assert_features(required_features)

print(espressomd.features())

# Interaction parameters (repulsive Lennard-Jones)
#############################################################

species = ["anion", "cation"]
types = {"anion": 0, "cation": 1}
charges = {"anion": -1.0, "cation": 1.0}
lj_sigmas = {"anion": 1.0, "cation": 1.0}
lj_epsilons = {"anion": 1.0, "cation": 1.0}
WCA_cut = 2.**(1. / 6.)
lj_cuts = {"anion": WCA_cut * lj_sigmas["anion"],
    "cation": WCA_cut * lj_sigmas["cation"]}
masses = {"anion": 1.0, "cation": 1.0}

# System parameters
#############################################################

try:
    nproc = int(os.environ.get("OMPI_COMM_WORLD_SIZE", 1))
    n_part_per_core = int(sys.argv[1])
    n_part = nproc * n_part_per_core
    bjerrum_length = float(sys.argv[2])
    mode = "benchmark"
    if len(sys.argv) == 4:
        assert sys.argv[3] == "--visualize"
        mode = "visualization"
        from espressomd import visualization
        from threading import Thread
except (ValueError, IndexError) as err:
    print(err.message)
    print("\nUsage: [mpiexec -np <N cores>] pypresso lj.py "
          "<N particles per core> <bjerrum_length> [--visualize]\n")
    exit(1)

measurement_steps = int(np.round(5e5 / n_part_per_core, -1))
assert(measurement_steps >= 50), \
    "{} steps per tick are too short".format(measurement_steps)
# volume of N spheres with radius r: N * (4/3*pi*r^3)
density = 0.25
lj_sig = (lj_sigmas["cation"] + lj_sigmas["anion"]) / 2
box_l = (n_part * 4. / 3. * np.pi * (lj_sig / 2.)**3 / density)**(1. / 3.)

# System
#############################################################
system = espressomd.System(box_l=3 * (box_l,))
system.cell_system.set_domain_decomposition(use_verlet_lists=True)
# PRNG seeds
#############################################################
system.random_number_generator_state = list(range(
    nproc * (system._get_PRNG_state_size() + 1)))
np.random.seed(1)
# Integration parameters
#############################################################
system.time_step = 0.01
system.cell_system.skin = 1.2
system.thermostat.turn_off()


#############################################################
#  Setup System                                             #
#############################################################

# Interaction setup
#############################################################

for i in range(len(species)):
    ion1 = species[i]
    for j in range(i, len(species)):
        ion2 = species[j]
        lj_sig = (lj_sigmas[ion1] + lj_sigmas[ion2]) / 2
        lj_cut = (lj_cuts[ion1] + lj_cuts[ion2]) / 2
        lj_eps = (lj_epsilons[ion1] * lj_epsilons[ion2])**(1. / 2.)
        system.non_bonded_inter[types[ion1],
            types[ion2]].lennard_jones.set_params(
                epsilon=lj_eps, sigma=lj_sig, cutoff=lj_cut, shift="auto")

# Particle setup
#############################################################

for i in range(0, n_part, len(species)):
    for t in species:
        system.part.add(pos=np.random.random(3) * system.box_l,
                        q=charges[t], type=types[t], mass=masses[t])

#############################################################
#  Warmup Integration                                       #
#############################################################

energy = system.analysis.energy()
print("Before Minimization: E_total = {}".format(energy['total']))
system.minimize_energy.init(f_max=1000, gamma=30.0,
                            max_steps=1000, max_displacement=0.05)
system.minimize_energy.minimize()
system.minimize_energy.minimize()
energy = system.analysis.energy()
print("After Minimization: E_total = {}".format(energy['total']))

print("Tune p3m")
p3m = electrostatics.P3M(prefactor=bjerrum_length, accuracy=1e-4)
system.actors.add(p3m)

system.thermostat.set_langevin(kT=1.0, gamma=1.0)

if mode == 'benchmark':
    report_path = 'benchmarks.csv'
    if '${CMAKE_BINARY_DIR}'[0] != '$':  # CMake variable
        report_path = '${CMAKE_BINARY_DIR}/benchmarks.csv.part'

    # print initial energies
    energies = system.analysis.energy()
    print(energies)

    # time integration loop
    print("Timing every {} steps".format(measurement_steps))
    main_tick = time()
    all_t = []
    for i in range(30):
        tick = time()
        system.integrator.run(measurement_steps)
        tock = time()
        t = (tock - tick) / measurement_steps
        print('step {}, time = {:.2e}, verlet: {:.2f}'
              .format(i, t, system.cell_system.get_state()["verlet_reuse"]))
        all_t.append(t)
    main_tock = time()
    # average time
    all_t = np.array(all_t)
    avg = np.average(all_t)
    ci = 1.96 * np.std(all_t) / np.sqrt(len(all_t) - 1)
    print("average: {:.3e} +/- {:.3e} (95% C.I.)".format(avg, ci))

    # print final energies
    energies = system.analysis.energy()
    print(energies)

    # write report
    report = ('"{script}","{arguments}",{cores},"{mpi}",{mean:.3e},'
              '{ci:.3e},{n},{dur:.1f},{E1:.5e},{E2:.5e},{E3:.5e}\n'.format(
                  script=os.path.basename(sys.argv[0]),
                  arguments=" ".join(map(str, sys.argv[1:])),
                  cores=nproc, dur=main_tock - main_tick, n=measurement_steps,
                  mpi="OMPI_COMM_WORLD_SIZE" in os.environ, mean=avg, ci=ci,
                  E1=system.analysis.energy()["total"],
                  E2=system.analysis.energy()["coulomb"],
                  E3=system.analysis.energy()["non_bonded"]))
    if not os.path.isfile(report_path):
        report = ('"script","arguments","cores","MPI","mean","ci",'
                  '"steps_per_tick","duration","E1","E2","E3"\n' + report)
    with open(report_path, 'a') as f:
        f.write(report)
else:
    # use visualizer
    visualizer = visualization.openGLLive(system)
    visualizer.run(1)
