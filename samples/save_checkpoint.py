"""
This samples demonstrates how to checkpoint a simulation.
"""

import espressomd

required_features = ["ELECTROSTATICS", "LENNARD_JONES"]
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
system.set_random_state_PRNG()
#system.seed = system.cell_system.get_state()['n_nodes'] * [1234]
np.random.seed(seed=system.seed)
system.time_step = 0.01
system.cell_system.skin = 0.4

checkpoint.register("system")

# test for "system.thermostat"
system.thermostat.set_langevin(kT=1.0, gamma=1.0)

# test for "system.non_bonded_inter"
lj_eps = 1.0
lj_sig = 1.0
lj_cut = 1.12246
lj_cap = 20

system.non_bonded_inter[0, 0].lennard_jones.set_params(
    epsilon=lj_eps, sigma=lj_sig,
    cutoff=lj_cut, shift="auto")
system.force_cap = lj_cap

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
