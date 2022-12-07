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

"""
Benchmark Lattice-Boltzmann fluid + Lennard-Jones particles.
"""
import espressomd
import espressomd.lb
import benchmarks
import numpy as np
import argparse

parser = argparse.ArgumentParser(description="Benchmark LB simulations. "
                                 "Save the results to a CSV file.")
parser.add_argument("--particles_per_core", metavar="N", action="store",
                    type=int, default=125, required=False,
                    help="Number of particles per core")
parser.add_argument("--box_l", action="store",
                    type=int, default=32, required=False,
                    help="Number of particles per core")
parser.add_argument("--lb_sites_per_particle", metavar="N_LB", action="store",
                    type=float, default=28, required=False,
                    help="Number of LB sites per particle")
parser.add_argument("--volume_fraction", metavar="FRAC", action="store",
                    type=float, default=0.03, required=False,
                    help="Fraction of the simulation box volume occupied by "
                    "particles (range: [0.01-0.74], default: 0.50)")
parser.add_argument("--single_precision", action="store_true", required=False,
                    help="Using single-precision floating point accuracy")
parser.add_argument("--output", metavar="FILEPATH", action="store",
                    type=str, required=False, default="benchmarks.csv",
                    help="Output file (default: benchmarks.csv)")

args = parser.parse_args()

# process and check arguments
n_iterations = 30
assert args.volume_fraction > 0, "volume_fraction must be a positive number"
assert args.volume_fraction < np.pi / (3 * np.sqrt(2)), \
    "volume_fraction exceeds the physical limit of sphere packing (~0.74)"

required_features = ["LENNARD_JONES", "WALBERLA"]
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
box_l = args.box_l
lb_grid = box_l
agrid = 1.
measurement_steps = 80

# System
#############################################################
system.box_l = 3 * (box_l,)

# Integration parameters
#############################################################
system.time_step = 0.01
system.cell_system.skin = 0.5
system.thermostat.turn_off()

print(f"LB shape: [{lb_grid}, {lb_grid}, {lb_grid}]")
print(f"LB agrid: {agrid:.3f}")

lbf = espressomd.lb.LBFluidWalberla(agrid=agrid, density=1, viscosity=1,
                                    tau=system.time_step,
                                    single_precision=args.single_precision)
system.actors.add(lbf)


# time integration loop
timings = benchmarks.get_timings(system, measurement_steps, n_iterations)

# average time
avg, ci = benchmarks.get_average_time(timings)
print(f"average: {1000 * avg:.1f} +/- {1000 * ci:.1f} ms (95% C.I.)")

# write report
benchmarks.write_report(args.output, n_proc, timings, measurement_steps)
