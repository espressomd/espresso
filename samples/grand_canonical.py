#
# Copyright (C) 2013-2019 The ESPResSo project
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
Perform a grand canonical simulation of a system in contact
with a salt reservoir while maintaining a constant chemical potential.
"""
epilog = """
Takes two command line arguments as input: 1) the reservoir salt
concentration in units of 1/sigma^3 and 2) the excess chemical
potential of the reservoir in units of kT.

The excess chemical potential of the reservoir needs to be determined
prior to running the grand canonical simulation using the script called
widom_insertion.py which simulates a part of the reservoir at the
prescribed salt concentration. Be aware that the reservoir excess
chemical potential depends on all interactions in the reservoir system.
"""
import numpy as np
import argparse

import espressomd
from espressomd import reaction_ensemble
from espressomd import electrostatics

required_features = ["P3M", "EXTERNAL_FORCES", "WCA"]
espressomd.assert_features(required_features)

parser = argparse.ArgumentParser(epilog=__doc__ + epilog)
parser.add_argument('cs_bulk', type=float,
                    help="bulk salt concentration [1/sigma^3]")
parser.add_argument('excess_chemical_potential', type=float,
                    help="excess chemical potential [kT] "
                         "(obtained from Widom's insertion method)")
args = parser.parse_args()

# System parameters
#############################################################

cs_bulk = args.cs_bulk
excess_chemical_potential_pair = args.excess_chemical_potential
box_l = 50.0

# Integration parameters
#############################################################
system = espressomd.System(box_l=[box_l, box_l, box_l])
np.random.seed(seed=42)

system.time_step = 0.01
system.cell_system.skin = 0.4
temperature = 1.0


#############################################################
#  Setup System                                             #
#############################################################

# Particle setup
#############################################################
# type 0 = HA
# type 1 = A-
# type 2 = H+

for i in range(int(cs_bulk * box_l**3)):
    system.part.add(pos=np.random.random(3) * system.box_l, type=1, q=-1)
    system.part.add(pos=np.random.random(3) * system.box_l, type=2, q=1)

wca_eps = 1.0
wca_sig = 1.0
types = [0, 1, 2]
for type_1 in types:
    for type_2 in types:
        system.non_bonded_inter[type_1, type_2].wca.set_params(
            epsilon=wca_eps, sigma=wca_sig)

RE = reaction_ensemble.ReactionEnsemble(
    temperature=temperature, exclusion_radius=2.0, seed=3)
RE.add_reaction(
    gamma=cs_bulk**2 * np.exp(excess_chemical_potential_pair / temperature),
    reactant_types=[], reactant_coefficients=[], product_types=[1, 2],
    product_coefficients=[1, 1], default_charges={1: -1, 2: +1})
print(RE.get_status())
system.setup_type_map([0, 1, 2])


RE.reaction(10000)

p3m = electrostatics.P3M(prefactor=2.0, accuracy=1e-3)
system.actors.add(p3m)
p3m_params = p3m.get_params()
for key, value in p3m_params.items():
    print("{} = {}".format(key, value))


# Warmup
#############################################################
# warmup integration (steepest descent)
warm_steps = 20
warm_n_times = 20
min_dist = 0.9 * wca_sig

# minimize energy using min_dist as the convergence criterion
system.integrator.set_steepest_descent(f_max=0, gamma=1e-3,
                                       max_displacement=0.01)
i = 0
while system.analysis.min_dist() < min_dist and i < warm_n_times:
    print("minimization: {:+.2e}".format(system.analysis.energy()["total"]))
    system.integrator.run(warm_steps)
    i += 1

print("minimization: {:+.2e}".format(system.analysis.energy()["total"]))
print()
system.integrator.set_vv()

# activate thermostat
system.thermostat.set_langevin(kT=temperature, gamma=.5, seed=42)

# MC warmup
RE.reaction(1000)

n_int_cycles = 10000
n_int_steps = 600
num_As = []
deviation = None
for i in range(n_int_cycles):
    RE.reaction(10)
    system.integrator.run(steps=n_int_steps)
    num_As.append(system.number_of_particles(type=1))
    if i > 2 and i % 50 == 0:
        print("HA", system.number_of_particles(type=0), "A-",
              system.number_of_particles(type=1), "H+",
              system.number_of_particles(type=2))
        concentration_in_box = np.mean(num_As) / box_l**3
        deviation = (concentration_in_box - cs_bulk) / cs_bulk * 100
        print("average num A {:.1f} +/- {:.1f}, average concentration {:.8f}, "
              "deviation to target concentration {:.2f}%".format(
                  np.mean(num_As),
                  np.sqrt(np.var(num_As, ddof=1) / len(num_As)),
                  concentration_in_box, deviation))
