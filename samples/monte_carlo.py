#
# Copyright (C) 2023 The ESPResSo project
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
Rapid prototyping of new Monte Carlo methods in Python.

This sample provides a re-implementation of the core functionality
of :ref:`reaction methods <Reaction methods>` in Python,
with a focus on the :ref:`constant pH <Constant pH>` and
:ref:`reaction ensemble <Reaction Ensemble>` methods.
See :ref:`Writing new Monte Carlo methods` for more details.
The sample simulates the acid-based titration of polyelectrolyte chains.

The sample is designed to run with the :ref:`kernprof` profiler attached:

.. code-block:: bash

    pypresso --kernprof monte_carlo.py --mode=core
    pypresso --kernprof monte_carlo.py --mode=python

"""

import numpy as np
import itertools
import argparse
import math
import time
import tqdm
import os

import espressomd
import espressomd.polymer
import espressomd.electrostatics
import espressomd.reaction_methods

required_features = ["P3M", "WCA"]
espressomd.assert_features(required_features)

parser = argparse.ArgumentParser(
    prog=f"pypresso --kernprof {os.path.basename(__file__)}",
    epilog=__doc__.lstrip().split("\n", 1)[0])
parser.add_argument("--mode", choices=["core", "python"], default="python",
                    help="use C++ (core) or Python (python) implementation")
parser.add_argument("--method", choices=["cph", "re"], default="cph",
                    help="use constant pH (cph) or reaction ensemble (re)")
args = parser.parse_args()


if "line_profiler" not in dir():
    def profile(func):
        def wrapper(*args, **kwargs):
            return func(*args, **kwargs)
        return wrapper


def factorial_Ni0_by_factorial_Ni0_plus_nu_i(nu_i, N_i0):
    value = 1.
    if nu_i > 0:
        value /= math.factorial(N_i0 + nu_i) // math.factorial(N_i0)
    elif nu_i < 0:
        value *= math.factorial(N_i0) // math.factorial(N_i0 + nu_i)
    return value


class SingleReaction:

    def __init__(self, **kwargs):
        self.reactant_types = kwargs["reactant_types"]
        self.reactant_coefficients = kwargs["reactant_coefficients"]
        self.product_types = kwargs["product_types"]
        self.product_coefficients = kwargs["product_coefficients"]
        self.gamma = kwargs["gamma"]
        self.accepted_moves = 0
        self.trial_moves = 0
        self.accumulator_potential_energy_difference_exponential = []
        self.nu_bar = sum(self.product_coefficients) - \
            sum(self.reactant_coefficients)

    def get_acceptance_rate(self):
        return self.accepted_moves / self.trial_moves

    def make_backward_reaction(self):
        return SingleReaction(
            gamma=1. / self.gamma, reactant_types=self.product_types,
            reactant_coefficients=self.product_coefficients,
            product_types=self.reactant_types,
            product_coefficients=self.reactant_coefficients)


class ReactionAlgorithm:

    def __init__(self, **kwargs):
        self.system = kwargs["system"]
        self.kT = kwargs["kT"]
        self.non_interacting_type = 100
        self.reactions = []
        self.particle_inside_exclusion_range_touched = False
        self.default_charges = {}
        self.m_empty_p_ids_smaller_than_max_seen_particle = []
        self.rng = np.random.default_rng(seed=kwargs["seed"])
        self.exclusion = espressomd.reaction_methods.ExclusionRadius(**kwargs)

    def set_non_interacting_type(self, type):
        self.non_interacting_type = type

    def add_reaction(self, **kwargs):
        """
        Set up a reaction in the forward and backward directions.
        """
        self.default_charges.update(kwargs["default_charges"])
        forward_reaction = SingleReaction(
            reactant_types=kwargs["reactant_types"],
            reactant_coefficients=kwargs["reactant_coefficients"],
            product_types=kwargs["product_types"],
            product_coefficients=kwargs["product_coefficients"],
            gamma=kwargs["gamma"])
        backward_reaction = forward_reaction.make_backward_reaction()
        self.reactions.append(forward_reaction)
        self.reactions.append(backward_reaction)

    @profile
    def reaction(self, steps):
        """
        Perform reaction steps. Chemical reactions are selected at random.
        """
        E_pot = self.system.analysis.potential_energy()
        random = self.rng.choice(len(self.reactions), size=steps, replace=True)
        for i in random:
            E_pot = self.generic_oneway_reaction(self.reactions[i], E_pot)

    @profile
    def generic_oneway_reaction(self, reaction, E_pot_old):
        """
        Carry out a generic one-way chemical reaction of the type
        `A + B + ... --> C + D + ...` and return the new potential
        energy if the trial move is accepted.
        """
        reaction.trial_moves += 1
        self.particle_inside_exclusion_range_touched = False
        if not self.all_reactant_particles_exist(reaction):
            return E_pot_old

        old_particle_numbers = self.save_old_particle_numbers(reaction)
        p_properties_old_state = self.make_reaction_attempt(reaction)

        if self.particle_inside_exclusion_range_touched:  # reject trial move
            self.restore_system(p_properties_old_state)
            self.particle_inside_exclusion_range_touched = False
            return E_pot_old

        E_pot_new = self.system.analysis.potential_energy()
        E_pot_diff = E_pot_new - E_pot_old
        bf = self.calculate_acceptance_probability(
            reaction, E_pot_diff, old_particle_numbers)
        reaction.accumulator_potential_energy_difference_exponential.append(
            math.exp(-E_pot_diff / self.kT))
        if self.rng.uniform() < bf:  # accept trial move
            self.delete_hidden_particles(p_properties_old_state)
            reaction.accepted_moves += 1
            return E_pot_new
        else:  # reject trial move
            self.restore_system(p_properties_old_state)
            return E_pot_old

    @profile
    def make_reaction_attempt(self, reaction):
        """
        Carry out a chemical reaction and save the old system state.
        """
        minimum_number_of_types = min(len(reaction.reactant_types),
                                      len(reaction.product_types))
        maximum_number_of_types = max(len(reaction.reactant_types),
                                      len(reaction.product_types))
        p_properties_old_state = {"changed_particles": [],
                                  "created_particles": [],
                                  "hidden_particles": []}

        for index in range(minimum_number_of_types):
            r_type = reaction.reactant_types[index]
            p_type = reaction.product_types[index]
            r_charge = self.default_charges[r_type]
            p_charge = self.default_charges[p_type]

            # change reactant particles to product particles
            size = min(reaction.reactant_coefficients[index],
                       reaction.product_coefficients[index])
            for random_pid in self.get_random_pids(r_type, size):
                p = self.system.part.by_id(random_pid)
                p.type = p_type
                p.q = p_charge
                p_properties_old_state["changed_particles"].append(
                    {"pid": random_pid, "type": r_type, "charge": r_charge})

            # measure stoichiometric excess
            delta_n = reaction.product_coefficients[index] - \
                reaction.reactant_coefficients[index]

            if delta_n > 0:
                # create product particles
                for _ in range(delta_n):
                    pid = self.create_particle(p_type)
                    self.check_exclusion_range(pid)
                    p_properties_old_state["created_particles"].append(
                        {"pid": pid, "type": p_type, "charge": p_charge})
            elif delta_n < 0:
                # hide reactant particles
                for random_pid in self.get_random_pids(r_type, -delta_n):
                    p_properties_old_state["hidden_particles"].append(
                        {"pid": random_pid, "type": r_type, "charge": r_charge})
                    self.check_exclusion_range(random_pid)
                    self.hide_particle(random_pid)

        # create/hide particles with non-corresponding replacement types
        for index in range(minimum_number_of_types, maximum_number_of_types):
            if len(reaction.product_types) < len(reaction.reactant_types):
                r_type = reaction.reactant_types[index]
                r_charge = self.default_charges[r_type]
                size = reaction.reactant_coefficients[index]
                # hide superfluous reactant particles
                for random_pid in self.get_random_pids(r_type, size):
                    p_properties_old_state["hidden_particles"].append(
                        {"pid": random_pid, "type": r_type, "charge": r_charge})
                    self.check_exclusion_range(random_pid)
                    self.hide_particle(random_pid)
            else:
                p_type = reaction.product_types[index]
                p_charge = self.default_charges[p_type]
                # create additional product particles
                for _ in range(reaction.product_coefficients[index]):
                    pid = self.create_particle(p_type)
                    self.check_exclusion_range(pid)
                    p_properties_old_state["created_particles"].append(
                        {"pid": pid, "type": p_type, "charge": p_charge})

        return p_properties_old_state

    def calculate_acceptance_probability(
            self, reaction, E_pot_diff, old_particle_numbers):
        raise NotImplementedError("Derived classes must implement this method")

    def check_exclusion_range(self, pid):
        self.particle_inside_exclusion_range_touched |= self.exclusion.check_exclusion_range(
            pid=pid)

    def all_reactant_particles_exist(self, reaction):
        for r_type in reaction.reactant_types:
            r_index = reaction.reactant_types.index(r_type)
            r_coef = reaction.reactant_coefficients[r_index]
            if self.system.number_of_particles(type=r_type) < r_coef:
                return False
        return True

    def get_random_pids(self, ptype, size):
        pids = self.system.analysis.call_method(
            "get_pids_of_type", ptype=ptype)
        indices = self.rng.choice(len(pids), size=size, replace=False)
        return [pids[i] for i in indices]

    def save_old_particle_numbers(self, reaction):
        old_particle_numbers = {}
        for r_type in reaction.reactant_types + reaction.product_types:
            old_particle_numbers[r_type] = self.system.number_of_particles(
                type=r_type)
        return old_particle_numbers

    def delete_created_particles(self, p_properties_old_state):
        for particle_info in p_properties_old_state["created_particles"]:
            self.system.part.by_id(particle_info["pid"]).remove()

    def delete_hidden_particles(self, p_properties_old_state):
        for particle_info in p_properties_old_state["hidden_particles"]:
            self.system.part.by_id(particle_info["pid"]).remove()

    def restore_system(self, p_properties_old_state):
        # restore properties of changed and hidden particles
        for particle_info in p_properties_old_state["changed_particles"] + \
                p_properties_old_state["hidden_particles"]:
            p = self.system.part.by_id(particle_info["pid"])
            p.type = particle_info["type"]
            p.q = particle_info["charge"]
        # destroy created particles
        self.delete_created_particles(p_properties_old_state)

    def hide_particle(self, pid):
        p = self.system.part.by_id(pid)
        p.type = self.non_interacting_type
        p.q = 0.

    def create_particle(self, ptype):
        if len(self.m_empty_p_ids_smaller_than_max_seen_particle) == 0:
            pid = self.system.part.highest_particle_id + 1
        else:
            pid = min(self.m_empty_p_ids_smaller_than_max_seen_particle)
        self.system.part.add(id=pid, type=ptype, q=self.default_charges[ptype],
                             pos=self.rng.random((3,)) * self.system.box_l,
                             v=self.rng.normal(size=3) * math.sqrt(self.kT))
        return pid


class ReactionEnsemble(ReactionAlgorithm):

    def calculate_acceptance_probability(
            self, reaction, E_pot_diff, old_particle_numbers):
        """
        Calculate the acceptance probability of a Monte Carlo move.
        """

        volume = self.system.volume()
        expr = math.exp(-E_pot_diff / self.kT)
        expr *= volume**reaction.nu_bar * reaction.gamma

        # factorial contribution of reactants
        for i in range(len(reaction.reactant_types)):
            nu_i = -reaction.reactant_coefficients[i]
            N_i0 = old_particle_numbers[reaction.reactant_types[i]]
            expr *= factorial_Ni0_by_factorial_Ni0_plus_nu_i(nu_i, N_i0)

        # factorial contribution of products
        for i in range(len(reaction.product_types)):
            nu_i = reaction.product_coefficients[i]
            N_i0 = old_particle_numbers[reaction.product_types[i]]
            expr *= factorial_Ni0_by_factorial_Ni0_plus_nu_i(nu_i, N_i0)

        return expr


class ConstantpHEnsemble(ReactionAlgorithm):

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.constant_pH = kwargs["constant_pH"]

    def add_reaction(self, **kwargs):
        kwargs["reactant_coefficients"] = [1]
        kwargs["product_coefficients"] = [1, 1]
        super().add_reaction(**kwargs)

    def calculate_acceptance_probability(
            self, reaction, E_pot_diff, old_particle_numbers):
        """
        Calculate the acceptance probability of a Monte Carlo move.
        """

        ln_bf = E_pot_diff - reaction.nu_bar * self.kT * math.log(10.) * (
            self.constant_pH + reaction.nu_bar * math.log10(reaction.gamma))

        factorial_expr = math.exp(-ln_bf / self.kT)

        # factorial contribution of reactants
        nu_i = -reaction.reactant_coefficients[0]
        N_i0 = old_particle_numbers[reaction.reactant_types[0]]
        factorial_expr *= factorial_Ni0_by_factorial_Ni0_plus_nu_i(nu_i, N_i0)

        # factorial contribution of products
        nu_i = reaction.product_coefficients[0]
        N_i0 = old_particle_numbers[reaction.product_types[0]]
        factorial_expr *= factorial_Ni0_by_factorial_Ni0_plus_nu_i(nu_i, N_i0)

        return factorial_expr


# System parameters
#############################################################
box_l = 35
system = espressomd.System(box_l=[box_l] * 3)
SEED = 42
N_acid = 50
N_chain = 5
sigma = 1.
epsilon = 1.
system.time_step = 0.01
system.cell_system.skin = 1.0
N_steps_MD = 1000
N_steps_MC = 50

# Reaction parameters
#############################################################
pKa = 7.
pH = 7.25
kT = 1.
friction = 1.
types = {
    "HA": 0,
    "A-": 1,
    "H+": 2,
}
charges = {
    "HA": 0.,
    "A-": -1.,
    "H+": +1.,
}
params = {
    "kT": kT,
    "exclusion_range": sigma,
    "seed": SEED,
}
if args.method == "cph":
    params["constant_pH"] = pH

# Setup
#############################################################
np.random.seed(seed=SEED)

positions = espressomd.polymer.linear_polymer_positions(
    n_polymers=N_chain, beads_per_chain=N_acid // N_chain,
    bond_length=sigma + 0.1, seed=SEED)

bond = espressomd.interactions.HarmonicBond(k=10., r_0=sigma + 0.1)
system.bonded_inter.add(bond)
for polymer_pos in positions:
    bond_partner = None
    for pos in polymer_pos:
        p = system.part.add(pos=pos, type=types["A-"], q=charges["A-"])
        if bond_partner:
            p.add_bond((bond, bond_partner))
        bond_partner = p

for _ in range(N_acid):
    pos = np.random.random(3) * system.box_l
    system.part.add(pos=pos, type=types["H+"], q=charges["H+"])

for type_pair in itertools.combinations_with_replacement(types.values(), 2):
    system.non_bonded_inter[type_pair[0], type_pair[1]].wca.set_params(
        epsilon=epsilon, sigma=sigma)

p3m = espressomd.electrostatics.P3M(
    prefactor=2., accuracy=1e-2, mesh=8, cao=3, verbose=False)
dh = espressomd.electrostatics.DH(
    prefactor=2., kappa=0., r_cut=0.2 * box_l)

# energy minimize the system
system.integrator.set_steepest_descent(
    f_max=0., gamma=0.1, max_displacement=0.1 * sigma)
system.integrator.run(200)
system.electrostatics.solver = p3m
system.integrator.run(1000)

# thermalize the system
system.integrator.set_vv()
system.thermostat.set_langevin(kT=kT, gamma=friction, seed=SEED)
system.integrator.run(1000)
system.electrostatics.solver = dh

if args.mode == "core":
    if args.method == "cph":
        ReactionMethod = espressomd.reaction_methods.ConstantpHEnsemble
    elif args.method == "re":
        ReactionMethod = espressomd.reaction_methods.ReactionEnsemble
elif args.mode == "python":
    if args.method == "cph":
        ReactionMethod = ConstantpHEnsemble
    elif args.method == "re":
        ReactionMethod = ReactionEnsemble
    params["system"] = system

# set up reaction method
RE = ReactionMethod(**params)
RE.set_non_interacting_type(type=max(types.values()) + 1)
if args.method == "cph":
    RE.add_reaction(
        gamma=10**-pKa,
        reactant_types=[types["HA"]],
        product_types=[types["A-"], types["H+"]],
        default_charges={types[name]: charges[name] for name in types.keys()})
elif args.method == "re":
    RE.add_reaction(
        gamma=1e-3,
        reactant_types=[types["HA"]],
        reactant_coefficients=[1],
        product_types=[types["A-"], types["H+"]],
        product_coefficients=[1, 1],
        default_charges={types[name]: charges[name] for name in types.keys()})
reaction = RE.reactions[0]
system.setup_type_map(type_list=list(types.values()))

# equilibrate the polyelectrolyte chains
for i in range(5):
    RE.reaction(steps=10 * N_steps_MC)
    system.integrator.run(N_steps_MD)


@profile
def sample_alpha(length):
    alpha_list = []
    for _ in tqdm.tqdm(range(length)):
        system.integrator.run(steps=N_steps_MD)
        RE.reaction(steps=N_steps_MC)
        alpha = system.number_of_particles(type=types["A-"]) / N_acid
        alpha_list.append(alpha)
    return alpha_list


sample_size = 100
tick = time.time()
alphas = sample_alpha(sample_size)
tock = time.time()

alpha_avg = np.mean(alphas)
alpha_err = 1.96 * np.sqrt(np.var(alphas) / len(alphas))
acceptance_rate = reaction.get_acceptance_rate()
print(f"acceptance rate = {100. * acceptance_rate:.0f}%")
print(f"alpha = {alpha_avg:.2f} +/- {alpha_err:.2f}")
print(f"runtime = {tock - tick:.2f}s")
