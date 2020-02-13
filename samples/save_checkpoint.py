# Copyright (C) 2010-2019 The ESPResSo project
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
"""
Basic usage of the checkpointing feature. Show how to write the state of:

* custom user variables.
* non-bonded interactions.
* particles.
* P3M parameters.
* thermostat.
"""

import espressomd

required_features = ["P3M", "WCA"]
espressomd.assert_features(required_features)

from espressomd import electrostatics
from espressomd import checkpointing

import numpy as np
import signal

checkpoint = checkpointing.Checkpoint(checkpoint_id="mycheckpoint")

if not checkpoint.checkpoint_signals:
    checkpoint.register_signal(signal.SIGINT)

# test for user data
myvar = "some script variable"
checkpoint.register("myvar")
myvar = "updated value"  # demo of how the register function works

# test for "system"
box_l = 10.7437

system = espressomd.System(box_l=[box_l] * 3)
np.random.seed(seed=42)
system.time_step = 0.01
system.cell_system.skin = 0.4

checkpoint.register("system")

# test for "system.thermostat"
system.thermostat.set_langevin(kT=1.0, gamma=1.0, seed=42)

# test for "system.non_bonded_inter"
wca_eps = 1.0
wca_sig = 1.0
wca_cap = 20

system.non_bonded_inter[0, 0].wca.set_params(epsilon=wca_eps, sigma=wca_sig)
system.force_cap = wca_cap

# test for "system.part"
n_part = 10
for i in range(n_part):
    system.part.add(id=i, pos=np.random.random(3) * system.box_l)

# test for "p3m"
for i in range(n_part // 2 - 1):
    system.part[2 * i].q = -1.0
    system.part[2 * i + 1].q = 1.0
p3m = electrostatics.P3M(prefactor=1.0, accuracy=1e-2)
system.actors.add(p3m)

# let's also register the p3m reference for easy access later
checkpoint.register("p3m")

checkpoint.save()
