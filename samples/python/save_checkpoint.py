import espressomd
from espressomd import checkpointing

checkpoint = checkpointing.Checkpointing("mycheckpoint", "./checkpointing/")

myvar = "some script variable"

checkpoint.register("myvar")

box_l = 10.7437

system = espressomd.System()
system.time_step = 0.01
system.cell_system.skin = 0.4
system.box_l = [box_l, box_l, box_l]

checkpoint.register("system")


lj_eps = 1.0
lj_sig = 1.0
lj_cut = 1.12246
lj_cap = 20

system.non_bonded_inter[0, 0].lennard_jones.set_params(
    epsilon=lj_eps, sigma=lj_sig,
    cutoff=lj_cut, shift="auto")
system.non_bonded_inter.set_force_cap(lj_cap)

checkpoint.register("system.non_bonded_inter")


checkpoint.save()
