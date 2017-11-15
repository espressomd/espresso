import espressomd
from espressomd import electrostatics
from espressomd import checkpointing

import numpy
import signal

checkpoint = checkpointing.Checkpointing(checkpoint_id="mycheckpoint")

# test for user data
myvar = "some script variable"
checkpoint.register("myvar")
myvar = "updated value"  # demo of how the register function works

skin = 0.4
checkpoint.register("skin")


# test for "system"
box_l = 10.7437

system = espressomd.System()
system.time_step = 0.01
system.cell_system.skin = skin
system.box_l = [box_l, box_l, box_l]

checkpoint.register("system")


# test for "system.thermostat"
system.thermostat.set_langevin(kT=1.0, gamma=1.0)

checkpoint.register("system.thermostat")


# test for "system.non_bonded_inter"
lj_eps = 1.0
lj_sig = 1.0
lj_cut = 1.12246
lj_cap = 20

system.non_bonded_inter[0, 0].lennard_jones.set_params(
    epsilon=lj_eps, sigma=lj_sig,
    cutoff=lj_cut, shift="auto")
system.force_cap = lj_cap

checkpoint.register("system.non_bonded_inter")


# test for "system.part"
n_part = 10
for i in range(n_part):
    system.part.add(id=i, pos=numpy.random.random(3) * system.box_l)

checkpoint.register("system.part")


# test for "p3m"
for i in range(n_part / 2 - 1):
    system.part[2 * i].q = -1.0
    system.part[2 * i + 1].q = 1.0
p3m = electrostatics.P3M(bjerrum_length=1.0, accuracy=1e-2)
system.actors.add(p3m)

checkpoint.register("p3m")

# signal.SIGINT: signal 2, is sent when ctrl+c is pressed
checkpoint.register_signal(signal.SIGINT)


checkpoint.save()
