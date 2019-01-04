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
from time import time
import argparse

parser = argparse.ArgumentParser(description="Benchmark P3M simulations. "
                                 "Save the results to a CSV file.")
parser.add_argument("--particles_per_core", metavar="N", action="store",
                    type=int, default=1000, required=False,
                    help="Number of particles in the simulation box")
parser.add_argument("--volume_fraction", metavar="FRAC", action="store",
                    type=float, default=0.25, required=False,
                    help="Fraction of the simulation box volume occupied by "
                    "particles (range: [0.01-0.74], default: 0.25)")
parser.add_argument("--bjerrum_length", metavar="LENGTH", action="store",
                    type=float, default=4., required=False,
                    help="Bjerrum length (default: 4)")
group = parser.add_mutually_exclusive_group()
group.add_argument("--output", metavar="FILEPATH", action="store",
                   type=str, required=False, default="benchmarks.csv",
                   help="Output file (default: benchmarks.csv)")
group.add_argument("--visualizer", action="store_true",
                   help="Starts the visualizer (for debugging purposes)")

args = parser.parse_args()

# process and check arguments
n_proc = int(os.environ.get("OMPI_COMM_WORLD_SIZE", 1))
n_part = n_proc * args.particles_per_core
measurement_steps = int(np.round(5e5 / args.particles_per_core, -1))
assert args.bjerrum_length > 0, "bjerrum_length must be a positive number"
assert args.volume_fraction > 0, "volume_fraction must be a positive number"
assert args.volume_fraction < np.pi / (3 * np.sqrt(2)), \
    "volume_fraction exceeds the physical limit of sphere packing (~0.74)"
if not args.visualizer:
    assert(measurement_steps >= 50), \
        "{} steps per tick are too short".format(measurement_steps)


import espressomd
from espressomd import thermostat
from espressomd import electrostatics
if args.visualizer:
    from espressomd import visualization
    from threading import Thread

required_features = ["ELECTROSTATICS", "LENNARD_JONES", "MASS"]
espressomd.assert_features(required_features)

print(espressomd.features())

# Interaction parameters (Lennard-Jones, Coulomb)
#############################################################

species = ["anion", "cation"]
types = {"anion": 0, "cation": 0}
charges = {"anion": -1.0, "cation": 1.0}
lj_sigmas = {"anion": 1.0, "cation": 1.0}
lj_epsilons = {"anion": 1.0, "cation": 1.0}
WCA_cut = 2.**(1. / 6.)
lj_cuts = {"anion": WCA_cut * lj_sigmas["anion"],
           "cation": WCA_cut * lj_sigmas["cation"]}
masses = {"anion": 1.0, "cation": 1.0}

# System parameters
#############################################################

# volume of N spheres with radius r: N * (4/3*pi*r^3)
lj_sig = (lj_sigmas["cation"] + lj_sigmas["anion"]) / 2
box_l = (n_part * 4. / 3. * np.pi * (lj_sig / 2.)**3
         / args.volume_fraction)**(1. / 3.)

# System
#############################################################
system = espressomd.System(box_l=3 * (box_l,))
system.cell_system.set_domain_decomposition(use_verlet_lists=True)
# PRNG seeds
#############################################################
system.random_number_generator_state = list(range(
    n_proc * (system._get_PRNG_state_size() + 1)))
# Integration parameters
#############################################################
system.time_step = 0.01
system.cell_system.skin = .4 
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
print("Before Minimization: E_total = {}".format(energy["total"]))
system.minimize_energy.init(f_max=1000, gamma=30.0,
                            max_steps=1000, max_displacement=0.05)
system.minimize_energy.minimize()
system.minimize_energy.minimize()
energy = system.analysis.energy()
print("After Minimization: E_total = {}".format(energy["total"]))


system.integrator.set_vv()
system.thermostat.set_langevin(kT=1.0, gamma=1.0)

system.integrator.run(min(3 * measurement_steps, 1000))
print("Tune skin: {}".format(system.cell_system.tune_skin(
    min_skin=0.4, max_skin=1.6, tol=0.05, int_steps=100)))
system.integrator.run(min(3 * measurement_steps, 3000))
print("Tune p3m")
p3m = electrostatics.P3M(prefactor=args.bjerrum_length, accuracy=1e-4)
system.actors.add(p3m)
system.integrator.run(min(3 * measurement_steps, 3000))
print("Tune skin: {}".format(system.cell_system.tune_skin(
    min_skin=1.0, max_skin=1.6, tol=0.05, int_steps=100)))


if not args.visualizer:
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
        print("step {}, time = {:.2e}, verlet: {:.2f}"
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
    cmd = " ".join(x for x in sys.argv[1:] if not x.startswith("--output"))
    report = ('"{script}","{arguments}",{cores},"{mpi}",{mean:.3e},'
              '{ci:.3e},{n},{dur:.1f},{E1:.5e},{E2:.5e},{E3:.5e}\n'.format(
                  script=os.path.basename(sys.argv[0]), arguments=cmd,
                  cores=n_proc, dur=main_tock - main_tick, n=measurement_steps,
                  mpi="OMPI_COMM_WORLD_SIZE" in os.environ, mean=avg, ci=ci,
                  E1=system.analysis.energy()["total"],
                  E2=system.analysis.energy()["coulomb"],
                  E3=system.analysis.energy()["non_bonded"]))
    if not os.path.isfile(args.output):
        report = ('"script","arguments","cores","MPI","mean","ci",'
                  '"steps_per_tick","duration","E1","E2","E3"\n' + report)
    with open(args.output, "a") as f:
        f.write(report)
else:
    # use visualizer
    visualizer = visualization.openGLLive(system)
    visualizer.run(1)
