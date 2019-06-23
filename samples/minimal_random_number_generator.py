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
import espressomd
import numpy as np
import sys

# System parameters
#############################################################

system = espressomd.System(box_l=[1.0, 1.0, 1.0])

n_nodes = 1  # for MPI
system.seed = np.random.randint(low=1, high=2**31 - 1, size=n_nodes)
# if no seed is provided espresso generates a seed
print("seed ", system.seed)
rng_state_read1 = system.random_number_generator_state
print("random number generator state read 1", rng_state_read1)

rng_state = []
for i in range(len(rng_state_read1)):
    rng_state.append(i)
system.random_number_generator_state = rng_state
rng_state_read2 = system.random_number_generator_state
print("random number generator state read 2", rng_state_read2)

system.set_random_state_PRNG()
rng_state_read3 = system.random_number_generator_state
print("random number generator state read 3", rng_state_read3)
