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
import os
import sys
import numpy as np
from time import time, sleep
import argparse

parser = argparse.ArgumentParser(description="Benchmark ferrofluid simulations. "
                                 "Save the results to a CSV file.")
parser.add_argument("--particles_per_core", metavar="N", action="store",
                    type=int, default=1000, required=False,
                    help="Number of particles in the simulation box")
parser.add_argument("--volume_fraction", metavar="FRAC", action="store",
                    type=float, default=0.1, required=False,
                    help="Fraction of the simulation box volume occupied by "
                    "particles (range: [0.01-0.74], default: 0.50)")
parser.add_argument("--dipole_moment", metavar="FRAC", action="store",
                    type=float, default=2**0.5, required=False,
                    help="Fraction of the simulation box volume occupied by "
                    "particles (range: [0.01-0.74], default: 0.50)")
group = parser.add_mutually_exclusive_group()
group.add_argument("--output", metavar="FILEPATH", action="store",
                   type=str, required=False, default="benchmarks.csv",
                   help="Output file (default: benchmarks.csv)")
group.add_argument("--visualizer", action="store_true",
                   help="Starts the visualizer (for debugging purposes)")

args = parser.parse_args()

# process and check arguments
measurement_steps = int(np.round(5e5 / args.particles_per_core, -2))
n_iterations = 20
assert args.volume_fraction > 0, "volume_fraction must be a positive number"
assert args.volume_fraction < np.pi / (3 * np.sqrt(2)), \
    "volume_fraction exceeds the physical limit of sphere packing (~0.74)"
assert args.dipole_moment > 0
if not args.visualizer:
    assert(measurement_steps >= 100), \
        "{} steps per tick are too short".format(measurement_steps)


import espressomd
from espressomd import magnetostatics
if args.visualizer:
    from espressomd import visualization
    from threading import Thread

required_features = ["LENNARD_JONES"]
espressomd.assert_features(required_features)

print(espressomd.features())

# System
#############################################################
system = espressomd.System(box_l=[1, 1, 1])

# Interaction parameters (Lennard-Jones)
#############################################################

lj_eps = 1.0  # LJ epsilon
lj_sig = 1.0  # particle diameter
lj_cut = lj_sig * 2**(1. / 6.)  # cutoff distance

# System parameters
#############################################################

n_proc = system.cell_system.get_state()['n_nodes']
n_part = n_proc * args.particles_per_core
# volume of N spheres with radius r: N * (4/3*pi*r^3)
box_l = (n_part * 4. / 3. * np.pi * (lj_sig / 2.)**3
         / args.volume_fraction)**(1. / 3.)

# System
#############################################################
system.box_l = 3 * (box_l,)

# PRNG seeds
#############################################################
# np.random.seed(1)

# Integration parameters
#############################################################
system.time_step = 0.01
system.cell_system.skin = 0.5
system.thermostat.turn_off()


#############################################################
#  Setup System                                             #
#############################################################

# Interaction setup
#############################################################
system.non_bonded_inter[0, 0].lennard_jones.set_params(
    epsilon=lj_eps, sigma=lj_sig, cutoff=lj_cut, shift="auto")

print("LJ-parameters:")
print(system.non_bonded_inter[0, 0].lennard_jones.get_params())

# Particle setup
#############################################################

for i in range(n_part):
    system.part.add(
        id=i,
        pos=np.random.random(3) *
        system.box_l,
        rotation=(1, 1, 1),
        dipm=args.dipole_moment)

#############################################################
#  Warmup Integration                                       #
#############################################################

system.integrator.set_steepest_descent(
    f_max=0,
    gamma=0.001,
    max_displacement=0.01)

# warmup
while system.analysis.energy()["total"] > 0.1 * n_part:
    print("minimization: {:.1f}".format(system.analysis.energy()["total"]))
    system.integrator.run(20)
print("minimization: {:.1f}".format(system.analysis.energy()["total"]))
print()
system.integrator.set_vv()

system.thermostat.set_langevin(kT=1.0, gamma=1.0, seed=42)

system.integrator.run(min(5 * measurement_steps, 60000))

system.actors.add(magnetostatics.DipolarP3M(prefactor=1, accuracy=1E-4))
print("Tune skin: {}".format(system.cell_system.tune_skin(
    min_skin=0.2, max_skin=1, tol=0.05, int_steps=100)))

system.integrator.run(min(5 * measurement_steps, 2500))

print(system.non_bonded_inter[0, 0].lennard_jones)

if not args.visualizer:
    # print initial energies
    energies = system.analysis.energy()
    print(energies)

    # time integration loop
    print("Timing every {} steps".format(measurement_steps))
    main_tick = time()
    all_t = []
    for i in range(n_iterations):
        tick = time()
        system.integrator.run(measurement_steps)
        tock = time()
        t = (tock - tick) / measurement_steps
        print("step {}, time = {:.2e}, verlet: {:.2f}, energy: {:.2e}"
              .format(i, t, system.cell_system.get_state()["verlet_reuse"],
                      system.analysis.energy()["total"]))
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
    report = ('"{script}","{arguments}",{cores},{mean:.3e},'
              '{ci:.3e},{n},{dur:.1f}\n'.format(
                  script=os.path.basename(sys.argv[0]), arguments=cmd,
                  cores=n_proc, dur=main_tock - main_tick, n=measurement_steps,
                  mean=avg, ci=ci))
    if not os.path.isfile(args.output):
        report = ('"script","arguments","cores","mean","ci",'
                  '"nsteps","duration"\n' + report)
    with open(args.output, "a") as f:
        f.write(report)
else:
    # use visualizer
    visualizer = visualization.openGLLive(system)

    def main_thread():
        while True:
            system.integrator.run(1)
            visualizer.update()
            sleep(1 / 60.)  # limit frame rate to at most 60 FPS

    t = Thread(target=main_thread)
    t.daemon = True
    t.start()
    visualizer.start()
