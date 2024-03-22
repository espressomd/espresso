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
import espressomd.electrostatics
import benchmarks
import numpy as np
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
parser.add_argument("--prefactor", metavar="PREFACTOR", action="store",
                    type=float, default=1., required=False,
                    help="P3M prefactor (default: 4)")
parser.add_argument("--gpu", action=argparse.BooleanOptionalAction,
                    default=False, required=False, help="Use GPU implementation")
group = parser.add_mutually_exclusive_group()
group.add_argument("--output", metavar="FILEPATH", action="store",
                   type=str, required=False, default="benchmarks.csv",
                   help="Output file (default: benchmarks.csv)")
group.add_argument("--visualizer", action="store_true",
                   help="Starts the visualizer (for debugging purposes)")

args = parser.parse_args()

# process and check arguments
measurement_steps = int(np.round(5e5 / args.particles_per_core, -1))
n_iterations = 30
assert args.prefactor > 0, "prefactor must be a positive number"
assert args.volume_fraction > 0, "volume_fraction must be a positive number"
assert args.volume_fraction < np.pi / (3 * np.sqrt(2)), \
    "volume_fraction exceeds the physical limit of sphere packing (~0.74)"
if not args.visualizer:
    assert measurement_steps >= 50, \
        f"{measurement_steps} steps per tick are too short"

required_features = ["P3M", "LENNARD_JONES"]
if args.gpu:
    required_features.append("CUDA")
espressomd.assert_features(required_features)

# make simulation deterministic
np.random.seed(42)

# System
#############################################################
system = espressomd.System(box_l=[1, 1, 1])

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

# System parameters
#############################################################
n_proc = system.cell_system.get_state()["n_nodes"]
n_part = n_proc * args.particles_per_core
# volume of N spheres with radius r: N * (4/3*pi*r^3)
lj_sig = (lj_sigmas["cation"] + lj_sigmas["anion"]) / 2
box_l = (n_part * 4. / 3. * np.pi * (lj_sig / 2.)**3
         / args.volume_fraction)**(1. / 3.)

# System
#############################################################
system.box_l = 3 * (box_l,)
system.cell_system.set_regular_decomposition(use_verlet_lists=True)

# Integration parameters
#############################################################
system.time_step = 0.01
system.cell_system.skin = 0.5

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
pid = 0
for i in range(0, n_part, len(species)):
    for t in species:
        system.part.add(pos=np.random.random(3) * system.box_l,
                        id=pid, q=charges[t], type=types[t])
        pid += 1

#  Warmup Integration
#############################################################

# warmup
benchmarks.minimize(system, n_part / 2.)

system.integrator.set_vv()
system.thermostat.set_langevin(kT=1.0, gamma=1.0, seed=42)

p3m_class = espressomd.electrostatics.P3M
if args.gpu:
    p3m_class = espressomd.electrostatics.P3MGPU

# tuning and equilibration
min_skin = 0.2
max_skin = 1.6
p3m_params = {"prefactor": args.prefactor, "accuracy": 1e-3}
p3m = p3m_class(**p3m_params)
print("Quick equilibration")
system.integrator.run(min(3 * measurement_steps, 1000))
print("Tune skin: {:.3f}".format(system.cell_system.tune_skin(
    min_skin=min_skin, max_skin=max_skin, tol=0.05, int_steps=100,
    adjust_max_skin=True)))
print("Equilibration")
system.integrator.run(min(3 * measurement_steps, 3000))
print("Tune p3m")
system.electrostatics.solver = p3m
print("Equilibration")
system.integrator.run(min(3 * measurement_steps, 3000))
print("Tune skin: {:.3f}".format(system.cell_system.tune_skin(
    min_skin=min_skin, max_skin=max_skin, tol=0.05, int_steps=100,
    adjust_max_skin=True)))


if args.visualizer:
    import espressomd.visualization
    visualizer = espressomd.visualization.openGLLive(system)
    visualizer.run(1)


# time integration loop
timings = benchmarks.get_timings(system, measurement_steps, n_iterations)

# average time
avg, ci = benchmarks.get_average_time(timings)
print(f"average: {avg:.3e} +/- {ci:.3e} (95% C.I.)")

# write report
benchmarks.write_report(args.output, n_proc, timings, measurement_steps)
