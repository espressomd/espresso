"""
This sample simulates the reaction ensemble. It also illustrates how the constant pH method can be used.
"""

#
# Copyright (C) 2013-2018 The ESPResSo project
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

import espressomd
from espressomd import code_info
from espressomd import analyze
from espressomd import integrate
from espressomd import reaction_ensemble


# System parameters
#############################################################
box_l = 35.0
volume=box_l**3

# Integration parameters
#############################################################
system = espressomd.System(box_l=[box_l] * 3)
np.random.seed(seed=343)

system.time_step = 0.02
system.cell_system.skin = 0.4



# Particle setup
#############################################################
# 2A +3B <-> 4C + 1D +3 E
N0 = 50  # number of reactant pairs

type_A=0
type_B=1
type_C=2
type_D=3
type_E=4
types=[type_A, type_B, type_C, type_D, type_E]

nu_A=-2
nu_B=-3
nu_C=4
nu_D=1
nu_E=3

nu_s=[nu_A, nu_B, nu_C, nu_D, nu_E]
nu_bar=np.sum(nu_s)

K=0.001
c_ref_in_mol_per_l=1.0 #reference concentration for which the reaction constant is reportet, 1mol/l
conversion_factor_from_1_div_sigma_cubed_to_mol_per_l=37.1 #for 1sigma=3.55Angstrom
c_ref_in_1_div_sigma_cubed=c_ref_in_mol_per_l/conversion_factor_from_1_div_sigma_cubed_to_mol_per_l 
Gamma=K*(c_ref_in_1_div_sigma_cubed)**nu_bar

for i in range(N0): #setup N0 pairs of reactants
    for j in range(abs(nu_A)):
        system.part.add(pos=np.random.random(3) * system.box_l, type=type_A)
    for j in range(abs(nu_B)):
        system.part.add(pos=np.random.random(3) * system.box_l, type=type_B)

RE = reaction_ensemble.ReactionEnsemble(temperature=1, exclusion_radius=0, seed=4) #exclusion radius 0 for ideal gas


RE.add_reaction(gamma=Gamma, reactant_types=[type_A, type_B], reactant_coefficients=[abs(nu_A),abs(nu_B)], product_types=[type_C, type_D, type_E], product_coefficients=[nu_C, nu_D, nu_E], default_charges={type_A: 0, type_B: 0, type_C: 0, type_D:0, type_E:0})
print(RE.get_status())

numbers={type_A: [], type_B:[], type_C:[], type_D:[], type_E:[]}

#warmup
RE.reaction(200)

for i in range(100):
    RE.reaction(10)
    for _type in types:
        numbers[_type].append(system.number_of_particles(type=_type))
    if(i % 100 == 0):
        print("iteration ", i)
        for _type in types:
            print("type: ", _type, numbers[_type][-1])

concentrations={type_A: np.mean(numbers[type_A])/volume*conversion_factor_from_1_div_sigma_cubed_to_mol_per_l, 
                type_B:np.mean(numbers[type_B])/volume*conversion_factor_from_1_div_sigma_cubed_to_mol_per_l, 
                type_C:np.mean(numbers[type_C])/volume*conversion_factor_from_1_div_sigma_cubed_to_mol_per_l, 
                type_D:np.mean(numbers[type_D])/volume*conversion_factor_from_1_div_sigma_cubed_to_mol_per_l, 
                type_E:np.mean(numbers[type_E])/volume*conversion_factor_from_1_div_sigma_cubed_to_mol_per_l}
                
std_err_concentrations={type_A: np.std(numbers[type_A],ddof=1)/len(numbers[type_A])/volume*conversion_factor_from_1_div_sigma_cubed_to_mol_per_l, 
                type_B:np.std(numbers[type_B],ddof=1)/len(numbers[type_B])/volume*conversion_factor_from_1_div_sigma_cubed_to_mol_per_l, 
                type_C:np.std(numbers[type_C],ddof=1)/len(numbers[type_C])/volume*conversion_factor_from_1_div_sigma_cubed_to_mol_per_l, 
                type_D:np.std(numbers[type_D],ddof=1)/len(numbers[type_D])/volume*conversion_factor_from_1_div_sigma_cubed_to_mol_per_l, 
                type_E:np.std(numbers[type_E],ddof=1)/len(numbers[type_E])/volume*conversion_factor_from_1_div_sigma_cubed_to_mol_per_l}                

print("concentrations sampled with the reaction ensemble (in mol/l)")
print(concentrations)
#print(std_err_concentrations)

import sympy

from sympy import *
c_A, c_B, c_C, c_D, c_E = symbols('c_A, c_B, c_C, c_D, c_E')
#K=c_C**nu_C*c_D**nu_D*c_E**nu_E*c_A**nu_A*c_B**nu_B
eq1 = K-((c_C/c_ref_in_mol_per_l)**nu_C*(c_D/c_ref_in_mol_per_l)**nu_D*(c_E/c_ref_in_mol_per_l)**nu_E*(c_A/c_ref_in_mol_per_l)**nu_A*(c_B/c_ref_in_mol_per_l)**nu_B)
eq2 = N0-(1.0/abs(nu_A)*c_A/conversion_factor_from_1_div_sigma_cubed_to_mol_per_l*volume + 1.0/nu_D*c_D/conversion_factor_from_1_div_sigma_cubed_to_mol_per_l*volume)
eq3 = c_A/c_B-float(nu_A)/float(nu_B)
eq4 = c_C/c_D-float(nu_C)/float(nu_D)
eq5 = c_C/c_E-float(nu_C)/float(nu_E)
initial_guess=[concentrations[type_A], concentrations[type_B], concentrations[type_C], concentrations[type_D], concentrations[type_E]]

print("expected concentrations in mol/l:")
c_A, c_B, c_C, c_D, c_E=nsolve((eq1, eq2, eq3, eq4, eq5), (c_A, c_B, c_C, c_D, c_E), initial_guess)
concentrations_sympy={type_A: c_A, type_B: c_B, type_C: c_C, type_D: c_D, type_E: c_E}



def test_concentrations(types):
    for type_ in types:
        if (abs(concentrations[type_]-concentrations_sympy[type_])>1e-3):
            raise RuntimeError("wrong concentration for type ", type_)
        else:
            print("concentration correct for type ", type_)

test_concentrations(types)

print("check equations solution reaction ensemble, the next numbers should equal 0")
print(K-((concentrations[type_C]/c_ref_in_mol_per_l)**nu_C*(concentrations[type_D]/c_ref_in_mol_per_l)**nu_D*(concentrations[type_E]/c_ref_in_mol_per_l)**nu_E*(concentrations[type_A]/c_ref_in_mol_per_l)**nu_A*(concentrations[type_B]/c_ref_in_mol_per_l)**nu_B))
print(N0-(1.0/abs(nu_A)*concentrations[type_A]/conversion_factor_from_1_div_sigma_cubed_to_mol_per_l*volume +1.0/nu_D*concentrations[type_D]/conversion_factor_from_1_div_sigma_cubed_to_mol_per_l*volume))
print(concentrations[type_A]/concentrations[type_B]-float(nu_A)/float(nu_B))
print(concentrations[type_C]/concentrations[type_D]-float(nu_C)/float(nu_D))
print(concentrations[type_C]/concentrations[type_E]-float(nu_C)/float(nu_E))