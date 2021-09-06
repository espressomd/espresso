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
parser.add_argument("--lb_sites_per_particle", metavar="N_LB", action="store",
                    type=float, default=28, required=False,
                    help="Number of LB sites per particle")
parser.add_argument("--volume_fraction", metavar="FRAC", action="store",
                    type=float, default=0.03, required=False,
                    help="Fraction of the simulation box volume occupied by "
                    "particles (range: [0.01-0.74], default: 0.50)")
parser.add_argument("--output", metavar="FILEPATH", action="store",
                    type=str, required=False, default="benchmarks.csv",
                    help="Output file (default: benchmarks.csv)")

args = parser.parse_args()

# process and check arguments
n_iterations = 30
assert args.volume_fraction > 0, "volume_fraction must be a positive number"
assert args.volume_fraction < np.pi / (3 * np.sqrt(2)), \
    "volume_fraction exceeds the physical limit of sphere packing (~0.74)"

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
lb_grid = int(2 * round((n_part * args.lb_sites_per_particle)**(1. / 3) / 2.))
agrid = box_l / lb_grid
measurement_steps = int(max(120**3 / lb_grid**3, 50))

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

system.part.add(pos=np.random.random((n_part, 3)) * system.box_l)

#  Warmup Integration
#############################################################

# warmup
benchmarks.minimize(system, n_part / 10.)

system.integrator.set_vv()
system.thermostat.set_langevin(kT=1.0, gamma=1.0, seed=42)

# tuning and equilibration
min_skin = 0.2
max_skin = 1.0
print("Tune skin: {:.3f}".format(system.cell_system.tune_skin(
    min_skin=min_skin, max_skin=max_skin, tol=0.05, int_steps=100)))
print("Equilibration")
system.integrator.run(500)
print("Tune skin: {:.3f}".format(system.cell_system.tune_skin(
    min_skin=min_skin, max_skin=max_skin, tol=0.05, int_steps=100)))
print("Equilibration")
system.integrator.run(500)


system.thermostat.turn_off()
print(f"LB shape: [{lb_grid}, {lb_grid}, {lb_grid}]")
print(f"LB agrid: {agrid:.3f}")
if hasattr(espressomd.lb, "LBFluid"):
    LBClass = espressomd.lb.LBFluid
elif hasattr(espressomd.lb, "LBFluidWalberla"):
    LBClass = espressomd.lb.LBFluidWalberla
else: 
    raise Exception("LB not built in")

lbf = LBClass(agrid=agrid, dens=1, visc=1, tau=system.time_step, kT=1, seed=1)
system.actors.add(lbf)
system.thermostat.set_lb(gamma=10, LB_fluid=lbf, seed=2)


# time integration loop
timings = benchmarks.get_timings(system, measurement_steps, n_iterations)

# average time
avg, ci = benchmarks.get_average_time(timings)
print(f"average: {avg:.3e} +/- {ci:.3e} (95% C.I.)")

# write report
benchmarks.write_report(args.output, n_proc, timings, measurement_steps)
