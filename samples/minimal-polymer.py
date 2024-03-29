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
Set up a linear polymer.
"""
import espressomd
import espressomd.interactions
import espressomd.polymer
import espressomd.io.writer.vtf  # pylint: disable=import-error

espressomd.assert_features(["WCA"])

import numpy as np

# System parameters
#############################################################

system = espressomd.System(box_l=[100, 100, 100])
np.random.seed(seed=42)

system.time_step = 0.01
system.cell_system.skin = 0.4
system.cell_system.set_n_square(use_verlet_lists=False)
outfile = open('polymer.vtf', 'w')

system.non_bonded_inter[0, 0].wca.set_params(epsilon=1, sigma=1)

fene = espressomd.interactions.FeneBond(k=10, d_r_max=2)
system.bonded_inter.add(fene)


positions = espressomd.polymer.linear_polymer_positions(
    n_polymers=1, beads_per_chain=50, bond_length=1.0, seed=3210)
previous_part = None
for pos in positions[0]:
    part = system.part.add(pos=pos)
    if previous_part:
        part.add_bond((fene, previous_part))
    previous_part = part

espressomd.io.writer.vtf.writevsf(system, outfile)


#############################################################
#      Warmup                                               #
#############################################################

# minimize energy using min_dist as the convergence criterion
system.integrator.set_steepest_descent(f_max=0, gamma=1e-3,
                                       max_displacement=0.01)
while system.analysis.min_dist() < 0.95:
    print(f"minimization: {system.analysis.energy()['total']:+.2e}")
    system.integrator.run(20)

print(f"minimization: {system.analysis.energy()['total']:+.2e}")
print()
system.integrator.set_vv()

# activate thermostat
system.thermostat.set_langevin(kT=1.0, gamma=1.0, seed=42)


#############################################################
#      Integration                                          #
#############################################################

print("simulating...")
t_steps = 1000
for t in range(t_steps):
    print(f"step {t + 1} of {t_steps}", end='\r', flush=True)
    system.integrator.run(10)
    espressomd.io.writer.vtf.writevcf(system, outfile)
outfile.close()
print()
