from __future__ import print_function

import espressomd
from espressomd import thermostat
from espressomd import interactions
from espressomd import integrate
from espressomd.virtual_sites import VirtualSitesRelative
from espressomd.io.writer import h5md

import numpy as np

system = espressomd.System(box_l=[1.0, 1.0, 1.0])
system.set_random_state_PRNG()
system.seed = system.cell_system.get_state()['n_nodes'] * [1234]

sim_ID="rigid"
h5md_file_name = "rigid_body.h5"

system.time_step = 0.01
skin=10.0
system.cell_system.skin = skin
box_l=100
system.box_l = [box_l, box_l, box_l]
system.thermostat.set_langevin(kT=1.0, gamma=1.0)

### Particle types
type_centre = 0
type_A = 1


for i in xrange(2):
    for j in xrange(2):
        system.non_bonded_inter[i, j].lennard_jones.set_params(
            epsilon=1., sigma=1.,
            cutoff=2.**(1. / 6.), shift="auto")


#############################################################
print("** Placing particles")
#############################################################


# start at p_id=1, and reserve p_id=0 for the centre bead.
p_id=1

branch_len=5
x=box_l*0.5
y=box_l*0.5
z=box_l*0.5


for n in xrange(branch_len):
    system.part.add(id=p_id, pos=[x+(n+1), y, z], type=type_A)
    p_id+=1
for n in xrange(branch_len):
    system.part.add(id=p_id, pos=[x-(n+1), y, z], type=type_A)
    p_id+=1
for n in xrange(branch_len):
    system.part.add(id=p_id, pos=[x, y+(n+1), z], type=type_A)
    p_id+=1
for n in xrange(branch_len):
    system.part.add(id=p_id, pos=[x, y-(n+1), z], type=type_A)
    p_id+=1
for n in xrange(branch_len):
    system.part.add(id=p_id, pos=[x, y, z+(n+1)], type=type_A)
    p_id+=1
for n in xrange(branch_len):
    system.part.add(id=p_id, pos=[x, y, z-(n+1)], type=type_A)
    p_id+=1


#############################################################
# rigid stuff
#############################################################
system.virtual_sites= VirtualSitesRelative(have_velocity=True)

#here we calculate the center of mass position (com) and the moment of inertia (momI) of the object
com=np.zeros(3)
momI=0
for i in range(1, 1+branch_len*6, 1):
    com+=system.part[i].pos
com/=(branch_len*6.)

max_dist=0
for i in range(1, 1+branch_len*6, 1):
    momI+=(com-system.part[i].pos)**2
    dist=np.sum((com-system.part[i].pos)**2)
    if dist>max_dist: max_dist=dist
max_dist=np.sqrt(max_dist)

# only change this for multiple GPUS
#system.min_global_cut = max_dist
print("max dist is {}".format(max_dist))
print("moment of intertia is", momI)

#place center bead
system.part.add(id=0, pos=com, mass=branch_len*6, rinertia=momI, rotation=[1,1,1], type=type_centre)

for i in range(1, 1+branch_len*6, 1):
    system.part[i].virtual=1
    system.part[i].vs_auto_relate_to(0)

############################################################



h5_file = h5md.H5md(filename=h5md_file_name, write_pos=True, write_vel=False,
                     write_force=False, write_species=True, write_mass=False,
                     write_charge=True, write_ordered=True)

h5_file.write()
h5_file.flush()

for frame in xrange(200):
    system.integrator.run(100)
    h5_file.write()
    h5_file.flush()
h5_file.close()

print("**Simulation finished")

