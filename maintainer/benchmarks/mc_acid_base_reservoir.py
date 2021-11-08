#
# Copyright (C) 2021 The ESPResSo project
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

import numpy as np
import pint
import time
import benchmarks
import espressomd
import espressomd.electrostatics
import espressomd.reaction_ensemble
import setuptools
import argparse

parser = argparse.ArgumentParser(description="Benchmark MC simulations in the grand-reaction ensemble. "
                                 "Save the results to a CSV file.")
parser.add_argument("--particles_per_core", metavar="N", action="store",
                    type=int, default=125, required=False,
                    help="Number of particles per core")
parser.add_argument("--p3m", action="store_true", required=False,
                    help="Use P3M instead of DH")
group = parser.add_mutually_exclusive_group()
group.add_argument("--output", metavar="FILEPATH", action="store",
                   type=str, required=False, default="benchmarks.csv",
                   help="Output file (default: benchmarks.csv)")
group.add_argument("--script-tune", action="store_true",
                   help="Tune the benchmark script parameters")
args = parser.parse_args()

# process and check arguments
assert args.particles_per_core >= 100, "you need to use at least 100 particles per core to avoid finite-size effects in the simulation"
espressomd.assert_features(['WCA', 'ELECTROSTATICS'])
assert setuptools.version.pkg_resources.packaging.specifiers.SpecifierSet('>=0.10.1').contains(pint.__version__), \
    f'pint version {pint.__version__} is too old: several numpy operations can cast away the unit'


def calc_ideal_alpha(pH, pKa):
    # ionization degree alpha calculated from the Henderson-Hasselbalch
    # equation for an ideal system
    return 1. / (1 + 10**(pKa - pH))


def calc_donnan_coefficient(c_acid, I_res, charge=-1):
    return charge * c_acid / (2 * I_res) + \
        np.sqrt((charge * c_acid / (2 * I_res))**2 + 1)


system = espressomd.System(box_l=[1.0] * 3)
n_proc = system.cell_system.get_state()['n_nodes']
n_part = n_proc * args.particles_per_core

# Important constants
ureg = pint.UnitRegistry()
TEMPERATURE = 298 * ureg.K
WATER_PERMITTIVITY = 78.5  # at 298 K
KT = TEMPERATURE * ureg.k_B
PARTICLE_SIZE = 0.355 * ureg.nm
MOLE = 6.022e23 / ureg.mole
BJERRUM_LENGTH = (ureg.e**2 / (4 * ureg.pi * ureg.eps0 *
                               WATER_PERMITTIVITY * KT)).to('nm')
pKw = 14

# Grand-reaction ensemble simulation: acid-base reaction inside the system + exchange of ions with the reservoir
# at pK=pH=7 and pSalt=3, we can neglect H+ and OH- ions. Then we need to implement only two reactions:
# 1. HA <=> A- + Na+
# 2. (emptyset) <=> Na+ + Cl-
# The script was tuned so that N_ACID=120 yields about 500 particles per
# simulation box (average total number of all particles)
N_ACID = (120 * n_part) // 500
C_ACID = 1e-2 * ureg('mol/L')
pKa = 7
pH = 7
pSalt = 2
# additional parameters of the simulation
RUN_INTEGRATION = True

# integration time step in reduced units
dt = 0.01

TOTAL_NUM_MC_STEPS = int(1e5)
NUM_SAMPLES = 100
INTEGRATION_STEPS_PER_SAMPLE = 100
assert TOTAL_NUM_MC_STEPS % NUM_SAMPLES == 0, \
    f"Total number of MC steps must be divisible by total number of samples, got {TOTAL_NUM_MC_STEPS} and {NUM_SAMPLES}"
MC_STEPS_PER_SAMPLE = TOTAL_NUM_MC_STEPS // NUM_SAMPLES

# definitions of reduced units
ureg.define(f'reduced_energy = {TEMPERATURE} * boltzmann_constant')
ureg.define(f'reduced_length = {PARTICLE_SIZE}')
ureg.define(f'reduced_charge = 1*e')
# store the values of some parameters in dimensionless reduced units, that
# will be later passed to ESPResSo
KT_REDUCED = KT.to('reduced_energy').magnitude
BJERRUM_LENGTH_REDUCED = BJERRUM_LENGTH.to('reduced_length').magnitude
PARTICLE_SIZE_REDUCED = PARTICLE_SIZE.to('reduced_length').magnitude
C_REF = 1.0 * ureg('mol/l')  # reference concentration for the pKa value
C_REF_REDUCED = (C_REF * MOLE).to('reduced_length^-3')
print(f"KT = {KT.to('reduced_energy'):4.3f}")
print(f"PARTICLE_SIZE = {PARTICLE_SIZE.to('reduced_length'):4.3f}")
print(f"BJERRUM_LENGTH = {BJERRUM_LENGTH.to('reduced_length'):4.3f}")

BOX_V = (N_ACID / (ureg.avogadro_constant * C_ACID)).to("nm^3")
BOX_L = np.cbrt(BOX_V)
BOX_L_REDUCED = BOX_L.to('reduced_length').magnitude
print("C_ACID:", C_ACID.to('mol/L'))
C_ACID_REDUCED_CORRECT = N_ACID / BOX_V.to('reduced_length^3')
print("1/C_ACID_REDUCED_CORRECT: ", 1 / C_ACID_REDUCED_CORRECT)

c_H_res = 10**(-pH) * ureg('mol/l')
c_H_reduced = (c_H_res * MOLE).to('reduced_length^-3')
c_salt_res = 10**(-pSalt) * ureg('mol/l')
c_acid_reduced = (C_ACID * MOLE).to('reduced_length^-3')
c_salt_res_reduced = (c_salt_res * MOLE).to('reduced_length^-3')
c_Cl_reduced = c_salt_res_reduced + c_acid_reduced
c_Na_reduced = c_salt_res_reduced
# Initial concentration of salt ions should be the same as in the reservoir
N_SALT = int((c_salt_res * BOX_V.to('L') * MOLE).magnitude)
print("N_SALT:", N_SALT)

ideal_alpha = calc_ideal_alpha(pH=pH, pKa=pKa)
print("ideal_alpha:", ideal_alpha)
max_donnan_coefficient = calc_donnan_coefficient(
    C_ACID.to('mol/l'), c_salt_res.to('mol/l'))
print("max_donnan:", max_donnan_coefficient)
print("log10(max_donnan):", np.log10(max_donnan_coefficient))

# equilibrium constant for the insertion of NaCl ion pairs
K_NaCl_reduced = c_Na_reduced.magnitude * c_Cl_reduced.magnitude

# equilibrium constants for the acid-base reaction
K_A = 10**(-pKa)
K_ANa_reduced = K_A * C_REF_REDUCED.magnitude * \
    (c_salt_res / c_H_res).magnitude

# #### Set the particle types and charges

# particle types of different species
TYPES = {
    "HA": 0,
    "A": 1,
    "Na": 2,
    "Cl": 3,
}
# particle charges of different species
CHARGES = {
    "HA": (0 * ureg.e).to("reduced_charge").magnitude,
    "A": (-1 * ureg.e).to("reduced_charge").magnitude,
    "Na": (+1 * ureg.e).to("reduced_charge").magnitude,
    "Cl": (-1 * ureg.e).to("reduced_charge").magnitude,
}

# ### Initialize the ESPResSo system
system.box_l = [BOX_L_REDUCED] * 3
system.time_step = dt
system.cell_system.skin = 2.0

# make simulation deterministic
np.random.seed(seed=10)

# ### Set up particles and bonded-interactions

# add the dissociated acid particles
system.part.add(pos=np.random.random((N_ACID, 3)) * BOX_L_REDUCED,
                type=[TYPES["A"]] * N_ACID,
                q=[CHARGES["A"]] * N_ACID)

# add the corresponding number of counter-ions (Na+)
system.part.add(pos=np.random.random((N_ACID, 3)) * BOX_L_REDUCED,
                type=[TYPES["Na"]] * N_ACID,
                q=[CHARGES["Na"]] * N_ACID)

# add salt ion pairs
system.part.add(pos=np.random.random((N_SALT, 3)) * BOX_L_REDUCED,
                type=[TYPES["Na"]] * N_SALT,
                q=[CHARGES["Na"]] * N_SALT)
system.part.add(pos=np.random.random((N_SALT, 3)) * BOX_L_REDUCED,
                type=[TYPES["Cl"]] * N_SALT,
                q=[CHARGES["Cl"]] * N_SALT)

print(f"The system was initialized with {len(system.part)} particles")

# ### Set up non-bonded-interactions

# set the WCA interaction betwee all particle pairs
for type_1 in TYPES.values():
    for type_2 in TYPES.values():
        if type_1 >= type_2:
            system.non_bonded_inter[type_1, type_2].wca.set_params(
                epsilon=1.0, sigma=1.0)

# relax the overlaps with steepest descent
system.integrator.set_steepest_descent(
    f_max=0, gamma=0.1, max_displacement=0.1)
system.integrator.run(20)
system.integrator.set_vv()  # to switch back to velocity Verlet

# add thermostat and short integration to let the system relax further
system.thermostat.set_langevin(kT=KT_REDUCED, gamma=1.0, seed=7)
system.integrator.run(steps=1000)

COULOMB_PREFACTOR = BJERRUM_LENGTH_REDUCED * KT_REDUCED
if args.p3m:
    coulomb = espressomd.electrostatics.P3M(prefactor=COULOMB_PREFACTOR,
                                            accuracy=1e-3)
else:
    KAPPA = np.sqrt(c_salt_res.to('mol/L').magnitude) / 0.304 / ureg.nm
    KAPPA_REDUCED = KAPPA.to('1/reduced_length').magnitude
    coulomb = espressomd.electrostatics.DH(prefactor=COULOMB_PREFACTOR,
                                           kappa=KAPPA_REDUCED,
                                           r_cut=1. / KAPPA_REDUCED)

system.actors.add(coulomb)

# ### Set up the constant pH ensemble using the reaction ensemble module
exclusion_radius = PARTICLE_SIZE_REDUCED
RE = espressomd.reaction_ensemble.ReactionEnsemble(
    kT=KT_REDUCED,
    exclusion_radius=exclusion_radius,
    seed=77
)
# this parameter helps speed up the calculation in an interacting system
RE.set_non_interacting_type(max(TYPES.values()) + 1)

RE.add_reaction(
    gamma=K_NaCl_reduced,
    reactant_types=[],
    reactant_coefficients=[],
    product_types=[TYPES["Na"], TYPES["Cl"]],
    product_coefficients=[1, 1],
    default_charges={TYPES["Na"]: CHARGES["Na"],
                     TYPES["Cl"]: CHARGES["Cl"]}
)

RE.add_reaction(
    gamma=K_ANa_reduced,
    reactant_types=[TYPES["HA"]],
    reactant_coefficients=[1],
    product_types=[TYPES["A"], TYPES["Na"]],
    product_coefficients=[1, 1],
    default_charges={TYPES["HA"]: CHARGES["HA"],
                     TYPES["A"]: CHARGES["A"],
                     TYPES["Na"]: CHARGES["Na"]}
)


def equilibrate_reaction(reaction_steps=1):
    RE.reaction(reaction_steps)


def report_progress(system, i, next_i):
    n_A = system.number_of_particles(type=TYPES['A'])
    n_Na = system.number_of_particles(type=TYPES['Na'])
    n_Cl = system.number_of_particles(type=TYPES['Cl'])
    n_All = len(system.part)
    if i == next_i:
        print(
            f"run {i:d} time {system.time:.3g} completed {i / NUM_SAMPLES * 100:.0f}%",
            f"instantaneous values: All {n_All:d}  Na {n_Na:d}  Cl {n_Cl:d}",
            f"A {n_A:d}  alpha {n_A / N_ACID:.3f}")
        if i == 0:
            next_i = 1
        else:
            next_i *= 2
    return next_i


# run a productive simulation and collect the data
start_time = time.time()
# pre-equilibrate the reaction to the new pH value
equilibrate_reaction(reaction_steps=10 * (N_ACID + 1))

# Main sampling loop
print(f"perform {NUM_SAMPLES:d} simulation runs")
if args.output:

    timings_MC = []
    timings_MD = []

    for i in range(NUM_SAMPLES):

        if MC_STEPS_PER_SAMPLE > 0:
            tick_MC = time.time()
            RE.reaction(MC_STEPS_PER_SAMPLE)
            tock_MC = time.time()

        t_MC = (tock_MC - tick_MC) / MC_STEPS_PER_SAMPLE
        timings_MC.append(t_MC)

        tick_MD = time.time()
        system.integrator.run(INTEGRATION_STEPS_PER_SAMPLE)
        tock_MD = time.time()

        t_MD = (tock_MD - tick_MD) / INTEGRATION_STEPS_PER_SAMPLE
        timings_MD.append(t_MD)

        energy = system.analysis.energy()["total"]
        verlet = system.cell_system.get_state()["verlet_reuse"]
        print(
            f"step {i}, time MD: {t_MD:.2e}, time MC: {t_MC:.2e}, verlet: {verlet:.2f}, energy: {energy:.2e}")

    # average time
    avg_MC, ci_MC = benchmarks.get_average_time(timings_MC)
    print(f"average MC: {avg_MC:.3e} +/- {ci_MC:.3e} (95% C.I.)")
    avg_MD, ci_MD = benchmarks.get_average_time(timings_MD)
    print(f"average MD: {avg_MD:.3e} +/- {ci_MD:.3e} (95% C.I.)")
    # write report
    benchmarks.write_report(args.output, n_proc, timings_MC, NUM_SAMPLES,
                            label='MC')
    benchmarks.write_report(args.output, n_proc, timings_MD, NUM_SAMPLES,
                            label='MD')
elif args.script_tune:
    next_i = 0
    n_As = []
    n_Alls = []

    # Use the piece of code below if you need to tune the script parameters
    for i in range(NUM_SAMPLES):
        if RUN_INTEGRATION:
            system.integrator.run(INTEGRATION_STEPS_PER_SAMPLE)
        RE.reaction(MC_STEPS_PER_SAMPLE)
        n_A = system.number_of_particles(type=TYPES['A'])
        n_As.append(n_A)
        n_All = len(system.part)
        n_Alls.append(n_All)
        next_i = report_progress(system, i, next_i)

    # calculate the averages at the end in order to check that we are having
    # the correct number of particles
    n_A_avg = np.average(np.array(n_As))
    n_All_avg = np.average(np.array(n_Alls))
    print(
        f"Average values: n_All {n_All_avg:.2f} , n_A {n_A_avg:.2f}, alpha {N_ACID:.4f}")
