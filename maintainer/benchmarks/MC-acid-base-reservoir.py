#!/usr/bin/env python
# coding: utf-8

# ionization degree alpha calculated from the Henderson-Hasselbalch equation for an ideal system
def calc_ideal_alpha(pH, pKa):
    return 1. / (1 + 10**(pKa - pH))

def calc_donnan_coefficient(c_acid, I_res, charge = -1):
    return charge*c_acid/(2*I_res) + np.sqrt( ( charge*c_acid/(2*I_res) )**2 + 1 )

import sys
import matplotlib.pyplot as plt
import numpy as np
import setuptools
import pint  # module for working with units and dimensions
import time
assert setuptools.version.pkg_resources.packaging.specifiers.SpecifierSet('>=0.10.1').contains(pint.__version__),   f'pint version {pint.__version__} is too old: several numpy operations can cast away the unit'

import espressomd
espressomd.assert_features(['WCA', 'ELECTROSTATICS'])
import espressomd.electrostatics
import espressomd.reaction_ensemble
import espressomd.polymer
#import espressomd.analyze
from espressomd.interactions import HarmonicBond

ureg = pint.UnitRegistry()

# Important constants
TEMPERATURE = 298 * ureg.K
WATER_PERMITTIVITY = 78.5 # at 298 K
KT = TEMPERATURE * ureg.k_B
PARTICLE_SIZE = 0.355 * ureg.nm
MOLE = 6.022e23 / ureg.mole
BJERRUM_LENGTH = (ureg.e**2 / (4 * ureg.pi * ureg.eps0 * WATER_PERMITTIVITY * KT)).to('nm')
pKw = 14
sys.exit()

# Input parameters of Curk, Yuan and Luijten
#N_ACID = 500 # the case studied in th
N_ACID = 100 # the test case compared with Espresso
BOX_L = 50*BJERRUM_LENGTH
pKa = 7
pH = 7
pSalt = 3
epsilon_reduced = 0.0 # value of epsilon in reduced units
# number of reaction samples per each pH value
t_max = 2.1e4 # 1e3 equilibration + 2e4 production
dt = 0.01

# additional parameters used by us
ideal = False # is set to true, all interactions are turned off to simualte an ideal system
time_per_sample = 10
steps_per_round = int(time_per_sample / dt)
NUM_SAMPLES = int(t_max / time_per_sample)
if(ideal):
    id_flag = 1
    # Simulate an interacting system with steric repulsion (Warning: it will be slower than without WCA!)
    USE_WCA = False
    # Simulate an interacting system with electrostatics 
    USE_ELECTROSTATICS = False
    # By default, we use the Debye-Huckel potential because it allows short runtime of the tutorial
    # You can also use the P3M method for full electrostatics (Warning: it will be much slower than DH!)
    USE_P3M = False
else:
    id_flag = 0
    USE_WCA = True
    USE_ELECTROSTATICS = True
    USE_P3M = True

if USE_ELECTROSTATICS:
    assert USE_WCA, "Use electrostatics only with a short-range repulsive potential. Otherwise, singularity occurs at zero distance."

sysname ="G-RxMC-id-{0:}_N-{1:.0f}_pK-{2:.1f}_pH-{3:.1f}_pS-{4:.1f}_e-{5:.1f}".format(
            id_flag, 
            N_ACID,
            pKa,
            pH,
            pSalt,
            epsilon_reduced,
            )
print("sysname:", sysname)
#sys.exit();

# definitions of reduced units
ureg.define(f'reduced_energy = {TEMPERATURE} * boltzmann_constant')
ureg.define(f'reduced_length = {PARTICLE_SIZE}')
ureg.define(f'reduced_charge = 1*e')
# store the values of some parameters in dimensionless reduced units, that will be later passed to ESPResSo
KT_REDUCED = KT.to('reduced_energy').magnitude
BJERRUM_LENGTH_REDUCED = BJERRUM_LENGTH.to('reduced_length').magnitude
PARTICLE_SIZE_REDUCED = PARTICLE_SIZE.to('reduced_length').magnitude
C_REF = 1.0* ureg('mol/l') # reference concentration for the pKa value
C_REF_REDUCED = ( C_REF*MOLE ).to('reduced_length^-3')
print("KT = {:4.3f}".format(KT.to('reduced_energy')))
print("PARTICLE_SIZE = {:4.3f}".format(PARTICLE_SIZE.to('reduced_length')))
print("BJERRUM_LENGTH = {:4.3f}".format(BJERRUM_LENGTH.to('reduced_length')))


#BOX_V = (N_ACID / (ureg.avogadro_constant * C_ACID)).to("nm^3")
#BOX_L = np.cbrt(BOX_V)
BOX_L_REDUCED = BOX_L.to('reduced_length').magnitude
BOX_V = BOX_L**3
C_ACID = N_ACID/MOLE/BOX_V
print( "C_ACID:", C_ACID.to('mol/L') )
C_ACID_REDUCED_CORRECT = N_ACID/BOX_V.to('reduced_length^3')
print("1/C_REDUCED_CORRECT: ", 1/C_ACID_REDUCED_CORRECT)

# Parameters that define properties of the system
# commented-out values are from teh testing version
#C_ACID = 0.1 * ureg('mol/l')
#N_ACID = 100 # desired number of salt ion pairs
#desired_alpha = 0.63 
#K_A = ( desired_alpha*C_ACID )**2 / ( (1-desired_alpha)*C_ACID  * C_REF )


c_ACID = 10**(-pH) * ureg('mol/l')
c_H = 10**(-pH) * ureg('mol/l')
c_H_reduced = (c_H * MOLE).to('reduced_length^-3')
c_OH = 10**(pH-pKw) * ureg ('mol/l')
c_OH_reduced = (c_OH * MOLE).to('reduced_length^-3')

c_neutralizer = np.abs(c_H - c_OH)
c_neutralizer_reduced = ( c_neutralizer * MOLE).to('reduced_length^-3')

C_SALT_RES = 10**(-pSalt) * ureg('mol/l')
C_SALT_RES_REDUCED = ( C_SALT_RES * MOLE ).to('reduced_length^-3')
if (c_H > c_OH):
    c_Cl_reduced = C_SALT_RES_REDUCED + c_neutralizer_reduced
    c_Na_reduced = C_SALT_RES_REDUCED 
    c_Cl = C_SALT_RES + c_neutralizer
    c_Na = C_SALT_RES 
else:
    c_Cl_reduced = C_SALT_RES_REDUCED 
    c_Na_reduced = C_SALT_RES_REDUCED + c_neutralizer_reduced
    c_Cl = C_SALT_RES 
    c_Na = C_SALT_RES + c_neutralizer

I_res1 = c_Na + c_H
I_res2 = c_Cl + c_OH
if ( np.abs(I_res1 - I_res2) / I_res1 < 1e-3):
    I_res = I_res1;
    print("I_res:", I_res)
else: 
    raise ValueError("Shoulbe be I_res1 = I_res2, got I_res[1,2]: {} {}\n".format(I_res1, I_res2) )

ideal_alpha = calc_ideal_alpha( pH=pH, pKa=pKa )
print("ideal_alpha:", ideal_alpha)
max_donnan_coefficient = calc_donnan_coefficient( C_ACID.to('mol/l'), I_res.to('mol/l') )
print("max_donnan:", max_donnan_coefficient)
print("log10(max_donnan):", np.log10(max_donnan_coefficient) )

# equilibrium constants for the insertion reactions
K_NaCl_reduced = c_Na_reduced.magnitude * c_Cl_reduced.magnitude 
K_NaOH_reduced = c_Na_reduced.magnitude * c_OH_reduced.magnitude 
K_HCl_reduced = c_H_reduced.magnitude * c_Cl_reduced.magnitude 
K_w_reduced = c_H_reduced.magnitude * c_OH_reduced.magnitude 

# equilibrium constants for the acid-base reactions
K_A = 10**(-pKa) 
K_A_reduced  = K_A * C_REF_REDUCED.magnitude 
K_AOH_reduced  = K_A_reduced / K_w_reduced
K_ANa_reduced  = K_A * C_REF_REDUCED.magnitude * (c_Na_reduced / c_H_reduced).magnitude
K_ACl_reduced  = K_A * C_REF_REDUCED.magnitude / (c_Cl_reduced * c_H_reduced).magnitude


#N_BLOCKS = 16  # number of block to be used in data analysis
#DESIRED_BLOCK_SIZE = 10  # desired number of samples per block

PROB_INTEGRATION = 1.0 # Probability of running MD integration after the reaction move. 
# This parameter changes the speed of convergence but not the limiting result
# to which the simulation converges


# #### Calculate the dependent parameters

#KAPPA = np.sqrt(C_SALT.to('mol/L').magnitude) / 0.304 / ureg.nm
#KAPPA_REDUCED = KAPPA.to('1/reduced_length').magnitude
#print(f"KAPPA = {KAPPA:.3f}")
#print(f"KAPPA_REDUCED = {KAPPA_REDUCED:.3f}")
#print(f"Debye_length: {1. / KAPPA:.2f} = {(1. / KAPPA).to('reduced_length'):.2f}")

print("Calculated values:")
print(f"BOX_L = {BOX_L:.3g} = {BOX_L.to('reduced_length'):.3g}")
print(f"BOX_V  = {BOX_V:.3g} = {BOX_V.to('reduced_length^3'):.3g}")
print(f"N_ACID = {N_ACID}")


# #### Set the particle types and charges

# particle types of different species
TYPES = {
    "HA": 0,
    "A": 1,
    "H": 2,
    "Na": 3,
    "OH": 4,
    "Cl": 5,
}
# particle charges of different species
CHARGES = {
    "HA": ( 0 * ureg.e).to("reduced_charge").magnitude,
    "A" : (-1 * ureg.e).to("reduced_charge").magnitude,
    "H" : (+1 * ureg.e).to("reduced_charge").magnitude,
    "OH": (-1 * ureg.e).to("reduced_charge").magnitude,
    "Na": (+1 * ureg.e).to("reduced_charge").magnitude,
    "Cl": (-1 * ureg.e).to("reduced_charge").magnitude,
}


# ### Initialize the ESPResSo system
system = espressomd.System(box_l=[BOX_L_REDUCED] * 3)
system.time_step = dt
system.cell_system.skin = 2.0
np.random.seed(seed=10)  # initialize the random number generator in numpy


# ### Set up particles and bonded-interactions
# we need to define bonds before creating polymers
hb = HarmonicBond(k=30, r_0=1.0)
system.bonded_inter.add(hb)

# create the polymer positions
polymers = espressomd.polymer.linear_polymer_positions(n_polymers=1,
                                                       beads_per_chain=N_ACID,
                                                       bond_length=0.9, 
                                                       seed=23)

# add the polymer particles composed of ionizable acid groups, initially in the ionized state
for polymer in polymers:
    prev_particle = None
    for position in polymer:
        p = system.part.add(pos=position, type=TYPES["A"], q=CHARGES["A"])
        if prev_particle:
            p.add_bond((hb, prev_particle))
        prev_particle = p

# add the corresponding number of H+ ions
system.part.add(pos=np.random.random((N_ACID, 3)) * BOX_L_REDUCED,
                type=[TYPES["H"]] * N_ACID,
                q=[CHARGES["H"]] * N_ACID)


## add salt ion pairs
#system.part.add(pos=np.random.random((N_SALT, 3)) * BOX_L_REDUCED,
#                type=[TYPES["Na"]] * N_SALT,
#                q=[CHARGES["Na"]] * N_SALT)
#system.part.add(pos=np.random.random((N_SALT, 3)) * BOX_L_REDUCED,
#                type=[TYPES["Cl"]] * N_SALT,
#                q=[CHARGES["Cl"]] * N_SALT)

print(f"The system was initialized with {len(system.part)} particles")


# ### Set up non-bonded-interactions
if USE_WCA:
    # set the WCA interaction betwee all particle pairs
    for type_1 in TYPES.values():
        for type_2 in TYPES.values():
            if type_1 >= type_2:
                system.non_bonded_inter[type_1, type_2].wca.set_params(epsilon=1.0, sigma=1.0)

    if(epsilon_reduced != 0.0):
        for type_1 in [TYPES['A'], TYPES['HA']]:
            for type_1 in [TYPES['A'], TYPES['HA']]:
                if type_1 >= type_2:
                    system.non_bonded_inter[type_1, type_2].lj.set_params(epsilon=epsilon, sigma=1.0, cut=2.5)
    # relax the overlaps with steepest descent
    system.integrator.set_steepest_descent(f_max=0, gamma=0.1, max_displacement=0.1)
    system.integrator.run(20)
    system.integrator.set_vv()  # to switch back to velocity Verlet

# add thermostat and short integration to let the system relax further
#system.thermostat.set_langevin(kT=KT_REDUCED, gamma=1.0, seed=7)
system.thermostat.set_langevin(kT=KT_REDUCED, gamma=1.0, seed=7)
system.integrator.run(steps=1000)

if USE_ELECTROSTATICS:
    COULOMB_PREFACTOR=BJERRUM_LENGTH_REDUCED * KT_REDUCED
    if USE_P3M:
        coulomb = espressomd.electrostatics.P3M(prefactor = COULOMB_PREFACTOR, 
                                                accuracy=1e-3)
    else:
        coulomb = espressomd.electrostatics.DH(prefactor = COULOMB_PREFACTOR, 
                                               kappa = KAPPA_REDUCED, 
                                               r_cut = 1. / KAPPA_REDUCED)
        
    system.actors.add(coulomb)
else:
    # this speeds up the simulation of dilute systems with small particle numbers
    system.cell_system.set_n_square()


# ### Set up the constant pH ensemble using the reaction ensemble module
exclusion_radius = PARTICLE_SIZE_REDUCED if USE_WCA else 0.0
RE = espressomd.reaction_ensemble.ReactionEnsemble(
    temperature=KT_REDUCED,
    exclusion_radius=exclusion_radius,
    seed=77
)
RE.set_non_interacting_type(len(TYPES)) # this parameter helps speed up the calculation in an interacting system

RE.add_reaction(
    gamma=K_NaCl_reduced,
    reactant_types=[],
    reactant_coefficients=[],
    product_types=[TYPES["Na"], TYPES["Cl"]],
    product_coefficients=[1,1],
    default_charges={TYPES["Na"]: CHARGES["Na"],
                     TYPES["Cl"]: CHARGES["Cl"],
                     }
)

RE.add_reaction(
    gamma=K_HCl_reduced,
    reactant_types=[],
    reactant_coefficients=[],
    product_types=[TYPES["H"], TYPES["Cl"]],
    product_coefficients=[1,1],
    default_charges={TYPES["H"]: CHARGES["H"],
                     TYPES["Cl"]: CHARGES["Cl"],
                     }
)

RE.add_reaction(
    gamma=K_NaOH_reduced,
    reactant_types=[],
    reactant_coefficients=[],
    product_types=[TYPES["Na"], TYPES["OH"]],
    product_coefficients=[1,1],
    default_charges={TYPES["Na"]: CHARGES["Na"],
                     TYPES["OH"]: CHARGES["OH"],
                     }
)

RE.add_reaction(
    gamma=K_w_reduced,
    reactant_types=[],
    reactant_coefficients=[],
    product_types=[TYPES["H"], TYPES["OH"]],
    product_coefficients=[1,1],
    default_charges={TYPES["H"]: CHARGES["H"],
                     TYPES["OH"]: CHARGES["OH"],
                     }
)

RE.add_reaction(
    gamma=K_A_reduced,
    reactant_types=[TYPES["HA"]],
    reactant_coefficients=[1],
    product_types=[TYPES["A"], TYPES["H"]],
    product_coefficients=[1,1],
    default_charges={TYPES["HA"]: CHARGES["HA"],
                     TYPES["A"]: CHARGES["A"],
                     TYPES["H"]: CHARGES["H"],
                     }
)

RE.add_reaction(
    gamma=K_AOH_reduced,
    reactant_types=[ TYPES["HA"], TYPES["OH"] ],
    reactant_coefficients=[1, 1],
    product_types=[ TYPES["A"] ],
    product_coefficients=[1],
    default_charges={TYPES["HA"]: CHARGES["HA"],
                     TYPES["A"]:  CHARGES["A"],
                     TYPES["OH"]: CHARGES["OH"],
                     }
)

RE.add_reaction(
    gamma=K_ANa_reduced,
    reactant_types=[TYPES["HA"]],
    reactant_coefficients=[1],
    product_types=[TYPES["A"], TYPES["Na"]],
    product_coefficients=[1,1],
    default_charges={TYPES["HA"]: CHARGES["HA"],
                     TYPES["A"]: CHARGES["A"],
                     TYPES["Na"]: CHARGES["Na"],
                     }
)

RE.add_reaction(
    gamma=K_ACl_reduced,
    reactant_types=[ TYPES["HA"], TYPES["Cl"] ],
    reactant_coefficients=[1, 1],
    product_types=[ TYPES["A"] ],
    product_coefficients=[1],
    default_charges={TYPES["HA"]: CHARGES["HA"],
                     TYPES["A"]:  CHARGES["A"],
                     TYPES["Cl"]: CHARGES["Cl"],
                     }
)

def equilibrate_reaction(reaction_steps=1):
    RE.reaction(reaction_steps)

def calc_observables(output_file, system, observables):
    values = []
    N_A = system.number_of_particles(type=TYPES["A"]) 
    N_HA = system.number_of_particles(type=TYPES["HA"]) 
    N_ACID = N_A + N_HA
    N_Na = system.number_of_particles(type=TYPES["Na"]) 
    N_Cl = system.number_of_particles(type=TYPES["Cl"]) 
    N_H = system.number_of_particles(type=TYPES["H"]) 
    N_OH = system.number_of_particles(type=TYPES["OH"]) 
    Rg_list =  system.analysis.calc_rg(chain_start=0, number_of_chains=1, chain_length = N_ACID)  
    Rg_value = Rg_list[0]*ureg('reduced_length')
    alpha = N_A / (N_ACID)
    for obs in observables:
        if obs == "time":
            output_file.write("{:.7g}  ".format(system.time))
        elif obs == "N_A":
            output_file.write("{:d}  ".format(N_A))
        elif obs == "N_HA":
            output_file.write("{:d}  ".format(N_HA))
        elif obs == "N_Na":
            output_file.write("{:d}  ".format(N_Na))
        elif obs == "N_Cl":
            output_file.write("{:d}  ".format(N_Cl))
        elif obs == "N_OH":
            output_file.write("{:d}  ".format(N_OH))
        elif obs == "N_H":
            output_file.write("{:d}  ".format(N_H))
        elif obs == "alpha":
            output_file.write("{:.3f}  ".format(alpha))
        elif obs == "BOX_V_dm3":
            output_file.write("{:.4g}  ".format( BOX_V.to('dm**3').magnitude) )
        elif obs == "Rg_nm":
            output_file.write("{:.4f}  ".format(Rg_value.to('nm').magnitude))
        else:
            raise ValueEror("Unknown observable " + obs + "\n")
    output_file.write("\n")
    output_file.flush()
    return N_A, alpha


def perform_sampling(output_file, observables, num_samples, reaction_steps, 
                     prob_integration=0.5, integration_steps=1000, ):
    alphas = [];
    obs_values = [];
    next_i = 0
    print("perform {:d} simulation runs".format(num_samples))
    for i in range(num_samples):
        if USE_WCA and np.random.random() < prob_integration:
            system.integrator.run(integration_steps)
        # we should do at least one reaction attempt per reactive particle
        RE.reaction(reaction_steps)
        n_A, alpha = calc_observables(
                output_file = output_file,
                system = system, 
                observables = observables, 
                )
#        current_n = system.number_of_particles(type=type_A) 
        alphas.append(alpha)
        if(i == next_i ):
            print("run {:d} time {:.3g} completed {:.2f}%  current_n_A {:d}  alpha {:.3f}".format(i, system.time, i/num_samples*100, n_A, alpha))
            if(i==0):
                i=1;
                next_i = 1;
            else:
                next_i = next_i * 2
    return np.array(alphas)

# run a productive simulation and collect the data
start_time = time.time()
equilibrate_reaction(reaction_steps=N_ACID + 1)  # pre-equilibrate the reaction to the new pH value

output_file = open(sysname + ".csv", "w")
observables =  ["time", "N_A", "N_HA", "N_Na", "N_Cl", "N_OH", "N_H", "alpha", "BOX_V_dm3", "Rg_nm"]
for obs in observables:
    output_file.write(obs + "  ")
output_file.write("\n")

# perform sampling
alphas = perform_sampling(output_file,
                 observables = observables,
                 num_samples=NUM_SAMPLES, 
                 reaction_steps=N_ACID + 1,
                 integration_steps = steps_per_round,
                 prob_integration=PROB_INTEGRATION)  # perform sampling / run production simulation
runtime = (time.time() - start_time) * ureg.s
mean_alpha = np.mean(alphas)
mean_n_A = N_ACID * mean_alpha
act_donnan_coefficient = calc_donnan_coefficient( mean_alpha*C_ACID.to('mol/l'), I_res.to('mol/l') )
ideal_alpha_with_donnan = calc_ideal_alpha(pH = pH+np.log10(act_donnan_coefficient), pKa = pKa)
print(f"runtime: {runtime:.2g}\n",
      f"measured number of A-: {mean_n_A:.2f}\n", 
      f"measured alpha: {mean_alpha:.2f}\n", 
      f"stdev: {np.std(alphas):.2f}\n", 
      f"ideal_alpha: {ideal_alpha}\n", 
      f"ideal_alpha with Donnan: {ideal_alpha_with_donnan}\n", 
      f"number of values: {alphas.shape}\n", 
      )
print("\nDONE")

sys.exit()

