#!/usr/bin/env python3
#
# Copyright (C) 2026 The ESPResSo project
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
Fixed-seed, serial, deterministic short trajectory for the ParticleStore
migration's bitwise-identity gate (spec section 5, phases 2-6).

Prints one line: '<mode> <sha256>' where the hash covers the bitwise
content of positions and forces after the run. Any change in this hash
between two commits of the same phase means the storage migration
altered numerics and is a bug (phases 2-6 are pure storage moves).
"""

import argparse
import hashlib

import numpy as np

import espressomd
import espressomd.electrostatics

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--mode", choices=["lj", "p3m"], required=True)
args = parser.parse_args()

np.random.seed(42)

system = espressomd.System(box_l=[12.0, 12.0, 12.0])
system.time_step = 0.001
system.cell_system.skin = 0.4

# 6x6x6 lattice with a small deterministic perturbation
positions = []
for i in range(6):
    for j in range(6):
        for k in range(6):
            positions.append([2.0 * i + 1.0, 2.0 * j + 1.0, 2.0 * k + 1.0])
positions = np.array(positions) + 0.1 * np.random.random((216, 3))

system.non_bonded_inter[0, 0].lennard_jones.set_params(
    epsilon=1.0, sigma=1.0, cutoff=2.5, shift="auto")

if args.mode == "lj":
    system.part.add(pos=positions)
else:
    charges = np.resize([1.0, -1.0], 216)
    system.part.add(pos=positions, q=charges)
    # fully pinned parameters: no tuning, fully deterministic. accuracy is a
    # required argument of the P3M script interface but is unused with
    # tune=False (mesh/cao/r_cut/alpha are set explicitly below).
    solver = espressomd.electrostatics.P3M(
        prefactor=1.0, accuracy=1e-4, mesh=[16, 16, 16], cao=6,
        r_cut=3.5, alpha=0.85, tune=False)
    system.electrostatics.solver = solver

system.integrator.run(20)

digest = hashlib.sha256()
digest.update(np.copy(system.part.all().pos).tobytes())
digest.update(np.copy(system.part.all().f).tobytes())
print(f"{args.mode} {digest.hexdigest()}")
