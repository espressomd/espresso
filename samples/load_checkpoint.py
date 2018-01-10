import espressomd
from espressomd import checkpointing

checkpoint = checkpointing.Checkpointing(checkpoint_id="mycheckpoint")
checkpoint.load()

# necessary for using e.g. system.actors, since explicit checkpointing of actors is not implemented yet
system = espressomd.System()

# test user variable
print "\n### user variable test ###"
print "myvar = {}".format(myvar)
print "skin = {}".format(skin)

# test "system"
print "\n### system test ###"
print "system.time = {}".format(system.time)
print "system.box_l = {}".format(system.box_l)
# system.cell_system not implemented yet, see sample script store_properties.py
system.cell_system.skin = skin

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

# test registered objects
# all objects that are registered when writing a checkpoint are automatically registered after loading this checkpoint
print "\n### checkpoint register test ###"
print "checkpoint.get_registered_objects() = {}".format(checkpoint.get_registered_objects())


# integrate system and finally save checkpoint
print "\n### Integrate until user presses ctrl+c ###"
print "Integrating..."
while True:
    system.integrator.run(1000)
