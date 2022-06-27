#
# Copyright (C) 2013-2022 The ESPResSo project
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

import espressomd
import espressomd.magnetostatics
import benchmarks
import numpy as np
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
                    help="Magnitude of the dipole moment (same for all particles)")
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
    assert measurement_steps >= 100, \
        f"{measurement_steps} steps per tick are too short"

required_features = ["LENNARD_JONES"]
espressomd.assert_features(required_features)

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

# Integration parameters
#############################################################
system.time_step = 0.01
system.cell_system.skin = 0.5
system.thermostat.turn_off()

# Interaction setup
#############################################################
system.non_bonded_inter[0, 0].lennard_jones.set_params(
    epsilon=lj_eps, sigma=lj_sig, cutoff=lj_cut, shift="auto")

# Particle setup
#############################################################
system.part.add(
    pos=np.random.random((n_part, 3)) * system.box_l,
    rotation=n_part * [(1, 1, 1)],
    dipm=n_part * [args.dipole_moment])

#  Warmup Integration
#############################################################

# warmup
benchmarks.minimize(system, n_part / 10.)

system.integrator.set_vv()
system.thermostat.set_langevin(kT=1.0, gamma=1.0, seed=42)

# tuning and equilibration
min_skin = 0.2
max_skin = 1.0
dp3m_params = {'prefactor': 1, 'accuracy': 1e-4}
print("Equilibration")
system.integrator.run(min(5 * measurement_steps, 60000))
dp3m = espressomd.magnetostatics.DipolarP3M(**dp3m_params)
system.actors.add(dp3m)
print("Tune skin: {:.3f}".format(system.cell_system.tune_skin(
    min_skin=min_skin, max_skin=max_skin, tol=0.05, int_steps=100)))
print("Equilibration")
system.integrator.run(min(5 * measurement_steps, 2500))


if args.visualizer:
    import time
    import threading
    import espressomd.visualization
    visualizer = espressomd.visualization.openGLLive(system)

    def main_thread():
        while True:
            system.integrator.run(1)
            visualizer.update()
            time.sleep(1 / 60.)  # limit frame rate to at most 60 FPS

    t = threading.Thread(target=main_thread)
    t.daemon = True
    t.start()
    visualizer.start()


# time integration loop
timings = benchmarks.get_timings(system, measurement_steps, n_iterations)

# average time
avg, ci = benchmarks.get_average_time(timings)
print(f"average: {avg:.3e} +/- {ci:.3e} (95% C.I.)")

# write report
benchmarks.write_report(args.output, n_proc, timings, measurement_steps)
