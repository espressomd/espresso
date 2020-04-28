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
Basic usage of the checkpointing feature. Show how to load the state of:

* custom user variables.
* non-bonded interactions.
* particles.
* P3M parameters.
* thermostat.
"""
# pylint: disable=undefined-variable
import espressomd

required_features = ["P3M", "WCA"]
espressomd.assert_features(required_features)

from espressomd import checkpointing
import numpy as np

checkpoint = checkpointing.Checkpoint(checkpoint_id="mycheckpoint")
checkpoint.load()

# print out actors

print("\n### current active actors ###")
for act in system.actors.active_actors:
    print(act)

# test user variable
print("\n### user variable test ###")
print("myvar = {}".format(myvar))

# test "system"
print("\n### system test ###")
print("system.time = {}".format(system.time))
print("system.box_l = {}".format(system.box_l))

# test "system.non_bonded_inter"
print("\n### system.non_bonded_inter test ###")
print("system.non_bonded_inter[0, 0].wca.get_params() = {}".format(
    system.non_bonded_inter[0, 0].wca.get_params()))

# test "system.part"
print("\n### system.part test ###")
print("system.part[:].pos = {}".format(system.part[:].pos))

# test "system.thermostat"
print("\n### system.thermostat test ###")
print("system.thermostat.get_state() = {}".format(
    system.thermostat.get_state()))

# test "p3m"
print("\n### p3m test ###")
print("p3m.get_params() = {}".format(p3m.get_params()))


# test registered objects
# all objects that are registered when writing a checkpoint are
# automatically registered after loading this checkpoint
print("\n### checkpoint register test ###")
print("checkpoint.get_registered_objects() = {}".format(
    checkpoint.get_registered_objects()))


# integrate system and finally save checkpoint
print("\n### Integrate until user presses ctrl+c ###")
print("Integrating...")

np.random.seed(seed=42)
while True:
    system.integrator.run(1000)
