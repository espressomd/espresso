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
"""
Guide for the constant pH ensemble. The modeled reaction is
:math:`2\\mathrm{A} + 3\\mathrm{B} \\leftrightarrow 4\\mathrm{C} + 1\\mathrm{H} + 3\\mathrm{E}`, where H represents a proton. In th constant pH ensemble an implicit coupling to a proton reservoir is assumed.
"""

import pprint
import numpy as np
import scipy.optimize

import espressomd
from espressomd import reaction_ensemble

# System parameters
#############################################################
box_l = 50.0
volume = box_l**3

# Integration parameters
#############################################################
system = espressomd.System(box_l=[box_l] * 3)
np.random.seed(seed=343)

system.time_step = 0.02
system.cell_system.skin = 0.4

# Particle setup
#############################################################
# 2A +3B <-> 4C + 1H +3 E
N0 = 50  # number of reactant pairs

type_A = 0
type_B = 1
type_C = 2
type_H = 3
type_E = 4
types = [type_A, type_B, type_C, type_H, type_E]
types_name = {type_A: 'A', type_B: 'B', type_C: 'C', type_H: 'H', type_E: 'E'}

nu_A = -2
nu_B = -3
nu_C = 4
nu_H = 5
nu_E = 3

nu_s = [nu_A, nu_B, nu_C, nu_H, nu_E]
nu_bar = np.sum(nu_s)

# reaction constant
K = 10**-3
# reference concentration for which the reaction constant is reported, 1 mol/l
c_ref_in_mol_per_l = 1.0
# simulation units: 1 sigma = 3.55 Angstrom
conversion_inv_sigma_cube_to_mol_per_l = 37.1
c_ref_in_1_div_sigma_cubed = c_ref_in_mol_per_l / \
    conversion_inv_sigma_cube_to_mol_per_l
pH = 1.0

# setup N0 pairs of reactants
for i in range(N0):
    for j in range(abs(nu_A)):
        system.part.add(pos=np.random.random(3) * system.box_l, type=type_A)
    for j in range(abs(nu_B)):
        system.part.add(pos=np.random.random(3) * system.box_l, type=type_B)

# use an exclusion radius of 0 to simulate an ideal gas
RE = reaction_ensemble.ConstantpHEnsemble(
    temperature=1, exclusion_radius=0, seed=4, conversion_factor_from_mol_per_l_to_1_div_sigma_cubed=1.0 / 37.1, pH=pH)

RE.add_reaction(
    gamma=K,
    reactant_types=[type_A, type_B],
    reactant_coefficients=[abs(nu_A), abs(nu_B)],
    product_types=[type_C, type_E, type_H],
    product_coefficients=[nu_C, nu_E, nu_H],
    default_charges={type_A: 0, type_B: 0, type_C: 0, type_H: 0, type_E: 0},
    product_index_protons=2)

pprint.pprint(RE.get_status())
RE.set_non_interacting_type(len(types))
numbers = {type_A: [], type_B: [], type_C: [], type_H: [], type_E: []}

# warmup
RE.reaction(300)

for i in range(100):
    RE.reaction(10)
    for _type in types:
        numbers[_type].append(system.number_of_particles(type=_type))

concentrations = {}
concentrations_95ci = {}
for ptype in types:
    concentrations[ptype] = np.mean(
        numbers[ptype]) / volume * conversion_inv_sigma_cube_to_mol_per_l
    concentrations_95ci[ptype] = 1.96 * np.std(numbers[ptype], ddof=1) / np.sqrt(len(
        numbers[ptype])) / volume * conversion_inv_sigma_cube_to_mol_per_l


def equations(variables):
    c_A, c_B, c_C, c_E = variables
    c_H = 10**(-pH) * c_ref_in_mol_per_l
    # K = c_C**nu_C * c_H**nu_H * c_E**nu_E * c_A**nu_A * c_B**nu_B
    eq1 = K - ((c_C / c_ref_in_mol_per_l)**nu_C * (c_H / c_ref_in_mol_per_l)**nu_H
               * (c_E / c_ref_in_mol_per_l)**nu_E * (c_A / c_ref_in_mol_per_l)**nu_A
               * (c_B / c_ref_in_mol_per_l)**nu_B)
    eq2 = N0 - (1.0 / abs(nu_A) * c_A / conversion_inv_sigma_cube_to_mol_per_l +
                1.0 / nu_C * c_C / conversion_inv_sigma_cube_to_mol_per_l) * volume
    eq3 = c_A / c_B - float(nu_A) / float(nu_B)
    eq4 = c_C / c_E - float(nu_C) / float(nu_E)
    return (eq1, eq2, eq3, eq4)


initial_guess = [
    concentrations[type_A],
    concentrations[type_B],
    concentrations[type_C],
    concentrations[type_E]]


c_A, c_B, c_C, c_E = scipy.optimize.fsolve(equations, initial_guess)

concentrations_numerical = {
    type_A: c_A,
    type_B: c_B,
    type_C: c_C,
    type_E: c_E}

concentrations_numerical[type_H] = 10**(-pH) * c_ref_in_mol_per_l

print("concentrations sampled with the reaction ensemble vs. analytical solutions:")
for ptype in types:
    if (ptype == type_H):
        print("Caution! Protons are only represented implicitly (i.e. not correctly if you count the protons in the simulation box)")
    print("  type {}: {:.4f} +/- {:.4f} mol/l (95% CI), expected: {:.10f} mol/l"
          .format(types_name[ptype],
                  concentrations[ptype],
                  concentrations_95ci[ptype],
                  concentrations_numerical[ptype]))

K_sim = ((concentrations[type_C] / c_ref_in_mol_per_l)**nu_C
         * (10**(-pH) * c_ref_in_mol_per_l / c_ref_in_mol_per_l)**nu_H
         * (concentrations[type_E] / c_ref_in_mol_per_l)**nu_E
         * (concentrations[type_A] / c_ref_in_mol_per_l)**nu_A
         * (concentrations[type_B] / c_ref_in_mol_per_l)**nu_B)
N0_sim = (1.0 / abs(nu_A) * concentrations[type_A] + 1.0 / nu_E *
          concentrations[type_E]) / conversion_inv_sigma_cube_to_mol_per_l * volume
print("properties of the simulated ensemble:")
print("  K_sim = {:.10f} mol/l, expected: {:.10f} mol/l".format(K_sim, K))
print("  diff in K = {:.4e} mol/l".format(K_sim - K))
print("  N0_sim = {:.1f}, expected: {}".format(N0_sim, N0))
