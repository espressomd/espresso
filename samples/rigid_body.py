from __future__ import print_function

import espressomd
from espressomd import thermostat
from espressomd import integrate
from espressomd.virtual_sites import VirtualSitesRelative
from espressomd.io.writer import h5md

import numpy as np

system = espressomd.System(box_l=[1.0, 1.0, 1.0])
system.set_random_state_PRNG()
system.seed = system.cell_system.get_state()['n_nodes'] * [1234]

sim_ID="rigid"
h5md_file_name = "rigid.h5"

system.time_step = 0.01
skin=10.0
system.cell_system.skin = skin
box_l=100
system.box_l = [box_l, box_l, box_l]
system.thermostat.set_langevin(kT=1.0, gamma=1.0)

### Particle types
type_centre = 0
type_A = 1


#############################################################
print("** Placing particles")
#############################################################


# start at p_id=1, and reserve p_id=0 for the centre bead.
p_id=1

branch_len=5
x=box_l*0.5
y=box_l*0.5
z=box_l*0.5

# place six branches, pointing +/-x +/-y and +/-z
for n in range(branch_len):
    system.part.add(pos=[x+(n+1), y, z], type=type_A, virtual=1)
for n in range(branch_len):
    system.part.add(pos=[x-(n+1), y, z], type=type_A, virtual=1)
for n in range(branch_len):
    system.part.add(pos=[x, y+(n+1), z], type=type_A, virtual=1)
for n in range(branch_len):
    system.part.add(pos=[x, y-(n+1), z], type=type_A, virtual=1)
for n in range(branch_len):
    system.part.add(pos=[x, y, z+(n+1)], type=type_A, virtual=1)
for n in range(branch_len):
    system.part.add(pos=[x, y, z-(n+1)], type=type_A, virtual=1)


#############################################################
# rigid stuff
#############################################################
system.virtual_sites= VirtualSitesRelative(have_velocity=True)

#here we calculate the center of mass position (com) and the moment of inertia (momI) of the object
com=np.zeros(3)
momI=0
for p in system.part:
    com+=p.pos
com/=(branch_len*6.)

max_dist=0
for p in system.part:
    momI+=(com-p.pos)**2
    dist=np.sum((p.pos)**2)
    if dist>max_dist: max_dist=dist
max_dist=np.sqrt(max_dist)

# only change this for multiple GPUS
#system.min_global_cut = max_dist
print("max dist is {}".format(max_dist))
print("moment of intertia is", momI)

#place center bead
p_center = system.part.add(pos=com, mass=branch_len*6, rinertia=momI, rotation=[1,1,1], type=type_centre)

# The virtual particles relate to the center one 
for p in system.part:
    if p.virtual:
        p.vs_auto_relate_to(p_center.id)

############################################################


h5_file = h5md.H5md(filename=h5md_file_name, write_pos=True, write_vel=False,
                     write_force=False, write_species=True, write_mass=False,
                     write_charge=True, write_ordered=True)

h5_file.write()
h5_file.flush()

for frame in range(200):
    system.integrator.run(100)
    h5_file.write()
    h5_file.flush()
h5_file.close()

print("**Simulation finished")

