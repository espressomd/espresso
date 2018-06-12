from __future__ import print_function

import espressomd
from espressomd import thermostat
from espressomd import integrate
from espressomd.virtual_sites import VirtualSitesRelative

import numpy as np

system = espressomd.System(box_l=[1.0, 1.0, 1.0])
system.set_random_state_PRNG()
system.seed = system.cell_system.get_state()['n_nodes'] * [1234]

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

system.virtual_sites= VirtualSitesRelative(have_velocity=True)

#here we calculate the center of mass position (com) and the moment of inertia (momI) of the object
com = system.analysis.center_of_mass(part_type=type_A)
print("center of mass is:", com)

# if using multiple nodes, we need to change min_global_cut to the largest deparation
if system.cell_system.get_state()['n_nodes'] >1:
    max_dist=0
    for p in system.part:
        dist=np.sum((p.pos-com)**2)
        if dist>max_dist: max_dist=dist
    max_dist=np.sqrt(max_dist)
    print("max dist is {}".format(max_dist))
    system.min_global_cut = max_dist

mat_I=system.analysis.moment_of_inertia_matrix(p_type=type_A)
#in this simple case, the cluster has principal axes aligned with the box
momI=[mat_I[0,0],mat_I[1,1], mat_I[2,2]]
print("moment of intertia is", momI)

#place center bead
p_center = system.part.add(pos=com, mass=branch_len*6, rinertia=momI, rotation=[1,1,1], type=type_centre)

# The virtual particles relate to the center one 
for p in system.part:
    if p.virtual:
        p.vs_auto_relate_to(p_center.id)





for frame in range(200):
    system.integrator.run(100)

print("**Simulation finished")

