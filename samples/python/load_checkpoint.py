import espressomd
from espressomd import checkpointing

checkpoint = checkpointing.Checkpointing(checkpoint_id="mycheckpoint")


checkpoint.load()

system = espressomd.System() #necessary for using e.g. system.actors, since explicit checkpointing of actors is not implemented yet

# test user variable
print "\n### user variable test ###"
print "myvar = {}".format(myvar)

# test "system"
print "\n### system test ###"
print "system.box_l = {}".format(system.box_l)
# system.cell_system not implemented yet, see sample script store_properties.py

# test "system.non_bonded_inter"
print "\n### system.non_bonded_inter test ###"
print "system.non_bonded_inter[0, 0].lennard_jones.get_params() = {}".format(system.non_bonded_inter[0, 0].lennard_jones.get_params())

# test "system.part"
print "\n### system.part test ###"
print "system.part[:].pos = {}".format(system.part[:].pos)

# test "system.thermostat"
print "\n### system.thermostat test ###"
print "system.thermostat.get_state() = {}".format(system.thermostat.get_state())

# test "p3m"
print "\n### p3m test ###"
print "p3m.get_params() = {}".format(p3m.get_params())
system.actors.add(p3m)
