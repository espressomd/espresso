#
# Copyright (C) 2013,2014 The ESPResSo project
#
# This file is part of ESPResSo.
#
# ESPResSo is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ESPResSo is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
from __future__ import print_function
import numpy as np
import sys

import espressomd
from espressomd import code_info
from espressomd import analyze
from espressomd import integrate
from espressomd.interactions import *
from espressomd import reaction_ensemble
from espressomd import electrostatics

# System parameters
#############################################################
cs_bulk=float(sys.argv[1])
excess_chemical_potential_pair=float(sys.argv[2]) #from widom insertion simulation of pair insertion
box_l=50.0

# Integration parameters
#############################################################
system = espressomd.System()
system.time_step = 0.01
system.cell_system.skin = 0.4
temperature=1.0
system.thermostat.set_langevin(kT=temperature, gamma=1.0)
system.cell_system.max_num_cells = 2744


#############################################################
#  Setup System                                             #
#############################################################

# Interaction setup
#############################################################
system.box_l = [box_l, box_l, box_l]

# Particle setup
#############################################################
# type 0 = HA
# type 1 = A-
# type 2 = H+

N0 = 50  # number of titratable units

for i in range(N0):
    system.part.add(id=i, pos=np.random.random(3) * system.box_l, type=1, q=-1)
for i in range(N0, 2 * N0):
    system.part.add(id=i, pos=np.random.random(3) * system.box_l, type=2, q=1)

lj_eps=1.0
lj_sig=1.0
lj_cut=2**(1.0/6)
types=[0,1,2]
for type_1 in types:
    for type_2 in types:
        system.non_bonded_inter[type_1, type_2].lennard_jones.set_params(
            epsilon=lj_eps, sigma=lj_sig,
            cutoff=lj_cut, shift="auto")

RE = reaction_ensemble.ReactionEnsemble(temperature=temperature, exclusion_radius=1.0, standard_pressure=1.0)
RE.add(equilibrium_constant=cs_bulk**2*np.exp(excess_chemical_potential_pair/temperature), reactant_types=[], reactant_coefficients=[], product_types=[1, 2], product_coefficients=[1, 1])
RE.set_default_charges(dictionary={"1": -1, "2": +1})
print(RE.get_status())
system.setup_type_map([0, 1, 2])


RE.reaction(10000) 

p3m = electrostatics.P3M(prefactor=0.9, accuracy=1e-3)
system.actors.add(p3m)
p3m_params = p3m.get_params()
for key in list(p3m_params.keys()):
    print("{} = {}".format(key, p3m_params[key]))


# Warmup
#############################################################
# warmup integration (with capped LJ potential)
warm_steps = 1000
warm_n_times = 20
# do the warmup until the particles have at least the distance min_dist
# set LJ cap
lj_cap = 20
system.force_cap=lj_cap

# Warmup Integration Loop
act_min_dist = system.analysis.min_dist()
i = 0
while (i < warm_n_times ):
    print(i, "warmup")
    RE.reaction(100)
    system.integrator.run(steps=warm_steps)
    i += 1
    #Increase LJ cap
    lj_cap = lj_cap + 10
    system.force_cap=lj_cap

# remove force capping
system.force_cap=0

#MC warmup
RE.reaction(1000) 

num_As=[]
for i in range(10000):
    RE.reaction(10)
    system.integrator.run(steps=300)
    num_As.append(system.number_of_particles(type=1))    
    if(i % 100 == 0):
        print("HA", system.number_of_particles(type=0), "A-",
              system.number_of_particles(type=1), "H+", system.number_of_particles(type=2))
        concentration_in_box=np.mean(num_As)/box_l**3
        print("average num A", np.mean(num_As), "+/-", np.sqrt(np.var(num_As,ddof=1)/len(num_As)), "average concentration", concentration_in_box , "deviation to target concentration", (concentration_in_box-cs_bulk)/cs_bulk*100, "%") 
        
        
              
