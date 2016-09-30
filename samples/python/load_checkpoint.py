from espressomd import checkpointing

checkpoint = checkpointing.Checkpointing("mycheckpoint", "./checkpointing/")


checkpoint.load()

print myvar
print system.box_l
print(system.non_bonded_inter[0, 0].lennard_jones.get_params())
