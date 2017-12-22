
#
# Copyright (C) 2013,2014,2015,2016 The ESPResSo project
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
import espressomd
from espressomd import thermostat
from espressomd import interactions
from espressomd import polymer
import numpy

# System parameters
#############################################################

system = espressomd.System()

# if no seed is provided espresso generates a seed

system.time_step = 0.01
system.cell_system.skin = 0.4
system.box_l = [100, 100, 100]
system.thermostat.set_langevin(kT=1.0, gamma=1.0)
system.cell_system.set_n_square(use_verlet_lists=False)

system.non_bonded_inter[0, 0].lennard_jones.set_params(
    epsilon=1, sigma=1,
    cutoff=2**(1. / 6), shift="auto")

fene = interactions.FeneBond(k=10, d_r_max=2)
system.bonded_inter.add(fene)

polymer.create_polymer(N_P=1, bond_length=1.0, MPC=50, bond=fene)

#############################################################
#      Integration                                          #
#############################################################

for i in range(20):
    system.integrator.run(1000)

    energies = system.analysis.energy()
    print(energies)
