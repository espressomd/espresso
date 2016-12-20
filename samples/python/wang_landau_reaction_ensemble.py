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
import espressomd._system as es
import espressomd
from espressomd import code_info
from espressomd import analyze
from espressomd import integrate
from espressomd.interactions import *
from espressomd import reaction_ensemble
import numpy as np

import sys
import os

print(code_info.features())
dev = "cpu"

# System parameters
#############################################################
box_l = 10

# Integration parameters
#############################################################
system = espressomd.System()
system.time_step = 0.5 #bigger timestep than usual is VERY important for hybrid monte carlo, at least 0.5 (or higher) for this system
system.cell_system.skin = 0.4
system.cell_system.max_num_cells = 2744


#############################################################
#  Setup System                                             #
#############################################################

# Interaction setup
#############################################################
system.box_l = [box_l, box_l, box_l]

# Particle setup
#############################################################
#type 0 = HA
#type 1 = A-
#type 2 = H+

N0 = 2 # number of titratable units
K_diss=0.0088

system.part.add(id=np.arange(N0) ,pos=np.random.random((N0,3))*box_l, type=np.zeros(N0))

h = HarmonicBond(r_0=0.0, k=10.0)
system.bonded_inter.add(h)
system.part[0].add_bond((h, 1))
bonded_energies=analyze.energy(system=system)["bonded"]

RE=reaction_ensemble.ReactionEnsemble(standard_pressure=0.00108, temperature=1, exclusion_radius=0)
RE.add(equilibrium_constant=K_diss,educt_types=[0],educt_coefficients=[1], product_types=[1,2], product_coefficients=[1,1])
RE.add(equilibrium_constant=1.0/K_diss, educt_types=[1,2], educt_coefficients=[1,1], product_types=[0],product_coefficients=[1])
RE.default_charges(dictionary={"0":0,"1":-1, "2":+1})
RE.print_status()

grand_canonical.setup([0,1,2])

## Wang-Landau stuff
RE.add_collective_variable_degree_of_association(associated_type=0,min=0, max=1, corresponding_acid_types=[0,1])
np.savetxt("temp_energy_file.dat",np.c_[[0,0.5,1],[0,0,0],[9,9,9]],header="nbar E_min E_max") # Note header may not be omitted since it is expected
RE.add_collective_variable_potential_energy(filename="temp_energy_file.dat",delta=0.5)
os.remove("temp_energy_file.dat")
RE.set_wang_landau_parameters(final_wang_landau_parameter=0.00001, wang_landau_steps=1, full_path_to_output_filename="WL.out", do_not_sample_reaction_partition_function=True, use_hybrid_monte_carlo=True)
RE.counter_ion_type=1



counter=1
while True:
	RE.do_reaction_wang_landau()
	
	##for preliminary run
#	counter+=1
#	RE.do_reaction_wang_landau()
#	RE.update_maximum_and_minimum_energies_at_current_state()
#	if counter%1000:
#		RE.write_out_preliminary_energy_run_results()
#	


#	_types=np.array(system.part[:].type)
#	print("HA", len(_types[_types==0]), "A-", len(_types[_types==1]), "H+", len(_types[_types==2]))
#	bonded_energies=analyze.energy(system=system)["bonded"]
#	print("E_bond", bonded_energies)
