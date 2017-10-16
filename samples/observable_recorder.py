import espressomd
import espressomd.observables

system = espressomd.System()
system.time_step = 0.01
system.cell_system.skin = 0

system.part.add(pos=[0, 0, 0], v=[1, 2, 3])
system.part.add(pos=[1, 0, 0], v=[2, 1, 3])
system.part.add(pos=[0, 1, 0], v=[3, 2, 1])

p = espressomd.observables.ComPosition(ids=[0, 1, 2])
p.auto_write_to(filename="test_observable_writer.dat")
system.auto_update_observables.add(p)

system.integrator.run(10)
