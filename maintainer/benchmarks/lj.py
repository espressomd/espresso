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

required_features = ["LENNARD_JONES"]
espressomd.assert_features(required_features)

print(espressomd.features())

# Interaction parameters (repulsive Lennard-Jones)
#############################################################

lj_eps = 1.0 # LJ epsilon
lj_sig = 1.0 # particle diameter
lj_cap = 20. # force cap
lj_cut = lj_sig * 2**(1. / 6.) # cutoff distance

# System parameters
#############################################################

try:
    nproc = int(os.environ.get("OMPI_COMM_WORLD_SIZE", 1))
    n_part_per_core = int(sys.argv[1])
    n_part = nproc * n_part_per_core
    if sys.argv[2] == 'gas':
        density = 0.02
    elif sys.argv[2] == 'liquid':
        density = 0.30
    else:
        density = float(sys.argv[2])
    assert density > 0, "density must be a positive number"
    assert density < np.pi / (3 * np.sqrt(2)), \
      "density exceeds the physical limit of sphere packing (~0.74)"
    mode = "benchmark"
    if len(sys.argv) == 4:
        assert sys.argv[3] == "--visualize"
        mode = "visualization"
        from espressomd import visualization
        from threading import Thread
except (ValueError, IndexError) as err:
    print(err.message)
    print("\nUsage: [mpiexec -np <N cores>] pypresso lj.py "
          "<N particles per core> <density> [--visualize]\n")
    exit(1)

measurement_steps = int(np.round(5e6 / n_part_per_core, -2))
assert(measurement_steps >= 100), \
  "{} steps per tick are too short".format(measurement_steps)
# volume of N spheres with radius r: N * (4/3*pi*r^3)
box_l = (n_part * 4. / 3. * np.pi * (lj_sig / 2.)**3 / density)**(1./3.)

# System
#############################################################
system = espressomd.System(box_l=3*(box_l,))
# PRNG seeds
#############################################################
system.random_number_generator_state = list(range(
      nproc * (system._get_PRNG_state_size() + 1)))
#system.random_number_generator_state = list(range(len(system.random_number_generator_state)))
np.random.seed(1)
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
system.force_cap = lj_cap

print("LJ-parameters:")
print(system.non_bonded_inter[0, 0].lennard_jones.get_params())

# Particle setup
#############################################################

for i in range(n_part):
    system.part.add(id=i, pos=np.random.random(3) * system.box_l)

#############################################################
#  Warmup Integration                                       #
#############################################################

system.integrator.set_steepest_descent(f_max=0, gamma=0.001, max_displacement=0.01)

# warmup
while system.analysis.energy()["total"] > 10 * n_part:
  print('minimization: {:.1f}'.format(system.analysis.energy()["total"]))
  system.integrator.run(10)
print()
system.integrator.set_vv()

system.thermostat.set_langevin(kT=1.0, gamma=1.0)

# tune skin
print('tune: {}'.format(system.cell_system.tune_skin(min_skin=0.2, max_skin=1, tol=0.05, int_steps=100)))
system.integrator.run(min(30 * measurement_steps, 60000))
print('tune: {}'.format(system.cell_system.tune_skin(min_skin=0.2, max_skin=1, tol=0.05, int_steps=100)))

print(system.non_bonded_inter[0, 0].lennard_jones)

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
        t = (tock-tick) / measurement_steps
        print('step {}, time = {:.2e}, verlet: {:.2f}'.format(i, t,
          system.cell_system.get_state()["verlet_reuse"]))
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
              E2=system.analysis.energy()["kinetic"],
              E3=system.analysis.energy()["non_bonded"]))
    if not os.path.isfile(report_path):
        report = ('"script","arguments","cores","MPI","mean","ci",'
                  '"steps_per_tick","duration","E1","E2","E3"\n' + report)
    with open(report_path, 'a') as f:
        f.write(report)
else:
    # use visualizer
    visualizer = visualization.openGLLive(system)

    def main_thread():
        while True:
            system.integrator.run(1)
            visualizer.update()
            sleep(1/60.) # limit framerate to at most 60 FPS

    t = Thread(target=main_thread)
    t.daemon = True
    t.start()
    visualizer.start()


