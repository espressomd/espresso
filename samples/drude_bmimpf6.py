# Copyright (C) 2010-2019 The ESPResSo project
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
"""
Particle polarization with cold Drude oscillators on a coarse-grained
simulation of the ionic liquid BMIM PF6.
"""

import time
import os
import numpy as np
import argparse

import espressomd
required_features = ["LENNARD_JONES", "P3M", "MASS", "ROTATION",
                     "ROTATIONAL_INERTIA", "VIRTUAL_SITES_RELATIVE",
                     "THOLE", "LANGEVIN_PER_PARTICLE"]
espressomd.assert_features(required_features)

from espressomd.minimize_energy import steepest_descent
from espressomd.electrostatics import P3M
from espressomd.interactions import ThermalizedBond, HarmonicBond
from espressomd import drude_helpers
from espressomd.virtual_sites import VirtualSitesRelative
import espressomd.visualization_opengl


print(__doc__ + """
The density and particle numbers are low for testing purposes. It writes the
xyz trajectory and RDF to the specified output path. Run with --help to see
available arguments (e.g. use --visual for visualization).
""")

parser = argparse.ArgumentParser(description='Drude LJ liquid')
parser.add_argument("--epsilon_r", nargs='?', default=1.0, type=float)
parser.add_argument("--mass_drude", nargs='?', default=0.8, type=float)
parser.add_argument("--walltime", nargs='?', default=1.0, type=float)
parser.add_argument("--drude", dest="drude", action="store_true")
parser.add_argument("--no-drude", dest="drude", action="store_false")
parser.add_argument("--thole", dest="thole", action="store_true")
parser.add_argument("--no-thole", dest="thole", action="store_false")
parser.add_argument("--intra_ex", dest="intra_ex", action="store_true")
parser.add_argument("--no-intra_ex", dest="intra_ex", action="store_false")
parser.add_argument("--visual", dest="visu", action="store_true")
parser.add_argument("--gpu", dest="gpup3m", action="store_true")
parser.add_argument("--path", nargs='?', default='./bmimpf6_bulk/', type=str)
parser.set_defaults(drude=True, thole=True,
                    intra_ex=True, visu=False, gpup3m=False)
args = parser.parse_args()

print("\nArguments:", args)

np.random.seed(42)
# NUM PARTICLES AND BOX
n_ionpairs = 100
n_part = n_ionpairs * 2
density_si = 0.5
# g/cm^3
rho_factor_bmim_pf6 = 0.003931
box_volume = n_ionpairs / rho_factor_bmim_pf6 / density_si
box_l = box_volume**(1. / 3.)
print("\n-->Ion pairs:", n_ionpairs, "Box size:", box_l)

system = espressomd.System(box_l=[box_l, box_l, box_l])
system.virtual_sites = VirtualSitesRelative()

if args.visu:
    d_scale = 0.988 * 0.5
    c_ani = [1, 0, 0, 1]
    c_dru = [0, 0, 1, 1]
    c_com = [0, 0, 0, 1]
    c_cat = [0, 1, 0, 1]
    visualizer = espressomd.visualization_opengl.openGLLive(
        system,
        background_color=[1, 1, 1],
        drag_enabled=True,
        ext_force_arrows=True,
        drag_force=10,
        draw_bonds=False,
        quality_particles=32,
        particle_coloring='type',
        particle_type_colors=[c_ani, c_cat, c_cat, c_cat,
                              c_com, c_dru, c_dru, c_dru, c_dru],
        particle_sizes=[
            0.5 * 5.06, 0.5 * 4.38, 0.5 * 3.41, 0.5 * 5.04, 0.1,
            d_scale * 5.06, d_scale * 4.38, d_scale * 3.41, d_scale * 5.04])

args.path = os.path.join(args.path, '')
if not os.path.exists(args.path):
    os.makedirs(args.path)

# TIMESTEP
fs_to_md_time = 1.0e-2
time_step_fs = 1.0
time_step_ns = time_step_fs * 1e-6
dt = time_step_fs * fs_to_md_time
system.time_step = dt

# TEMPERATURE
SI_temperature = 353.0
gamma_com = 1.0
kb_kjmol = 0.0083145
temperature_com = SI_temperature * kb_kjmol

# COULOMB PREFACTOR (elementary charge)^2 / (4*pi*epsilon_0) in Angstrom *
# kJ/mol
coulomb_prefactor = 1.67101e5 * kb_kjmol / args.epsilon_r

# DRUDE/TOTAL MASS
mass_tot = 100.0
mass_core = mass_tot - args.mass_drude
mass_red_drude = args.mass_drude * mass_core / mass_tot

# SPRING CONSTANT DRUDE
k_drude = 4184.0  # in kJ/mol/A^2
# Period of free oscillation: T_spring = 2Pi/w; w = sqrt(k_d/m_d)
T_spring = 2.0 * np.pi * np.sqrt(args.mass_drude / k_drude)

# TEMP DRUDE
SI_temperature_drude = 1.0
temperature_drude = SI_temperature_drude * kb_kjmol
gamma_drude = mass_red_drude / T_spring

# CELLSYSTEM
system.cell_system.skin = 0.4

# FORCEFIELD
types = {"PF6": 0, "BMIM_C1": 1, "BMIM_C2": 2, "BMIM_C3":
         3, "BMIM_COM": 4, "PF6_D": 5, "BMIM_C1_D": 6, "BMIM_C2_D": 7, "BMIM_C3_D": 8}
charges = {"PF6": -0.78, "BMIM_C1": 0.4374,
           "BMIM_C2": 0.1578, "BMIM_C3": 0.1848, "BMIM_COM": 0}
polarizations = {"PF6": 4.653, "BMIM_C1":
                 5.693, "BMIM_C2": 2.103, "BMIM_C3": 7.409}
masses = {"PF6": 144.96, "BMIM_C1": 67.07,
          "BMIM_C2": 15.04, "BMIM_C3": 57.12, "BMIM_COM": 0}
masses["BMIM_COM"] = masses["BMIM_C1"] + \
    masses["BMIM_C2"] + masses["BMIM_C3"]
lj_sigmas = {"PF6": 5.06, "BMIM_C1": 4.38,
             "BMIM_C2": 3.41, "BMIM_C3": 5.04, "BMIM_COM": 0}
lj_epsilons = {"PF6": 2.56, "BMIM_C1": 2.56,
               "BMIM_C2": 0.36, "BMIM_C3": 1.83, "BMIM_COM": 0}
shortTypes = ["PF6", "IM", "BU", "ME", "COM", "PF6_D", "IM_D", "BU_D", "ME_D"]
lj_types = ["PF6", "BMIM_C1", "BMIM_C2", "BMIM_C3"]

cutoff_sigmafactor = 2.5
lj_cuts = {}
for t in lj_sigmas:
    lj_cuts[t] = cutoff_sigmafactor * lj_sigmas[t]

system.min_global_cut = 3.5


def combination_rule_epsilon(rule, eps1, eps2):
    if rule == "Lorentz":
        return (eps1 * eps2)**0.5
    else:
        return ValueError("No combination rule defined")


def combination_rule_sigma(rule, sig1, sig2):
    if rule == "Berthelot":
        return (sig1 + sig2) * 0.5
    else:
        return ValueError("No combination rule defined")


# Lennard-Jones interactions parameters
for i in range(len(lj_types)):
    for j in range(i, len(lj_types)):
        s = [lj_types[i], lj_types[j]]
        lj_sig = combination_rule_sigma(
            "Berthelot", lj_sigmas[s[0]], lj_sigmas[s[1]])
        lj_cut = combination_rule_sigma(
            "Berthelot", lj_cuts[s[0]], lj_cuts[s[1]])
        lj_eps = combination_rule_epsilon(
            "Lorentz", lj_epsilons[s[0]], lj_epsilons[s[1]])

        system.non_bonded_inter[types[s[0]], types[s[1]]].lennard_jones.set_params(
            epsilon=lj_eps, sigma=lj_sig, cutoff=lj_cut, shift="auto")

# Place Particles
rid = 0
anion_ids = []
cation_sites_ids = []
cation_com_ids = []
cation_c1_ids = []
cation_c2_ids = []
cation_c3_ids = []

for i in range(n_ionpairs):
    system.part.add(id=rid, type=types["PF6"], pos=np.random.random(3) * box_l,
                    q=charges["PF6"], mass=masses["PF6"])
    anion_ids.append(rid)
    if args.drude:
        rid += 2
    else:
        rid += 1

for i in range(n_ionpairs):
    pos_com = np.random.random(3) * box_l
    system.part.add(
        id=rid, type=types["BMIM_COM"], pos=pos_com,
        mass=masses["BMIM_COM"], rinertia=[646.284, 585.158, 61.126],
        temp=0, gamma=0, rotation=[1, 1, 1])
    cation_com_ids.append(rid)
    com_id = rid
    rid += 1
    system.part.add(id=rid, type=types["BMIM_C1"],
                    pos=pos_com + [0, -0.527, 1.365], q=charges["BMIM_C1"])
    system.part[rid].vs_auto_relate_to(com_id)
    cation_c1_ids.append(rid)
    cation_sites_ids.append(rid)
    if args.drude:
        rid += 2
    else:
        rid += 1
    system.part.add(id=rid, type=types["BMIM_C2"],
                    pos=pos_com + [0, 1.641, 2.987], q=charges["BMIM_C2"])
    system.part[rid].vs_auto_relate_to(com_id)
    cation_c2_ids.append(rid)
    cation_sites_ids.append(rid)
    if args.drude:
        rid += 2
    else:
        rid += 1
    system.part.add(id=rid, type=types["BMIM_C3"],
                    pos=pos_com + [0, 0.187, -2.389], q=charges["BMIM_C3"])
    system.part[rid].vs_auto_relate_to(com_id)
    cation_c3_ids.append(rid)
    cation_sites_ids.append(rid)
    if args.drude:
        rid += 2
    else:
        rid += 1

# ENERGY MINIMIZATION
print("\n-->E minimization")
print("Before:", system.analysis.energy()["total"])
n_max_steps = 100000
steepest_descent(system, f_max=5.0, gamma=0.01, max_steps=n_max_steps,
                 max_displacement=0.01)
print("After:", system.analysis.energy()["total"])

# THERMOSTAT
if not args.drude:
    system.thermostat.set_langevin(
        kT=temperature_com,
        gamma=gamma_com,
        seed=42)

# ELECTROSTATICS
if args.gpup3m:
    from espressomd.electrostatics import P3MGPU
    print("\n-->Tune P3M GPU")
    p3m = P3MGPU(prefactor=coulomb_prefactor, accuracy=1e-3)
else:
    print("\n-->Tune P3M CPU")
    p3m = P3M(prefactor=coulomb_prefactor, accuracy=1e-3)

system.actors.add(p3m)

if args.drude:
    print("-->Adding Drude related bonds")
    thermalized_dist_bond = ThermalizedBond(
        temp_com=temperature_com, gamma_com=gamma_com,
        temp_distance=temperature_drude, gamma_distance=gamma_drude,
        r_cut=min(lj_sigmas.values()) * 0.5, seed=123)
    harmonic_bond = HarmonicBond(k=k_drude, r_0=0.0, r_cut=1.0)
    system.bonded_inter.add(thermalized_dist_bond)
    system.bonded_inter.add(harmonic_bond)

    for i in anion_ids:
        drude_helpers.add_drude_particle_to_core(
            system, harmonic_bond, thermalized_dist_bond, system.part[i],
            i + 1, types["PF6_D"], polarizations["PF6"],
            args.mass_drude, coulomb_prefactor)
    for i in cation_c1_ids:
        drude_helpers.add_drude_particle_to_core(
            system, harmonic_bond, thermalized_dist_bond, system.part[i],
            i + 1, types["BMIM_C1_D"], polarizations["BMIM_C1"],
            args.mass_drude, coulomb_prefactor)
    for i in cation_c2_ids:
        drude_helpers.add_drude_particle_to_core(
            system, harmonic_bond, thermalized_dist_bond, system.part[i],
            i + 1, types["BMIM_C2_D"], polarizations["BMIM_C2"],
            args.mass_drude, coulomb_prefactor)
    for i in cation_c3_ids:
        drude_helpers.add_drude_particle_to_core(
            system, harmonic_bond, thermalized_dist_bond, system.part[i],
            i + 1, types["BMIM_C3_D"], polarizations["BMIM_C3"],
            args.mass_drude, coulomb_prefactor)

    drude_helpers.setup_and_add_drude_exclusion_bonds(system)

    if args.thole:
        print("-->Adding Thole interactions")
        drude_helpers.add_all_thole(system)

    if args.intra_ex:
        # SETUP BONDS ONCE
        print("-->Adding intramolecular exclusions")
        drude_helpers.setup_intramol_exclusion_bonds(
            system,
            [types["BMIM_C1_D"], types["BMIM_C2_D"], types["BMIM_C3_D"]],
            [types["BMIM_C1"], types["BMIM_C2"], types["BMIM_C3"]],
            [charges["BMIM_C1"], charges["BMIM_C2"], charges["BMIM_C3"]])

        # ADD SR EX BONDS PER MOLECULE
        for i in cation_c1_ids:
            drude_helpers.add_intramol_exclusion_bonds(
                system, [i + 1, i + 3, i + 5], [i, i + 2, i + 4])

print("\n-->Short equilibration with smaller time step")
system.time_step = 0.1 * fs_to_md_time
system.integrator.run(1000)
system.time_step = time_step_fs * fs_to_md_time


print("\n-->Timing")
start = time.time()
n_timing_steps = 1000
system.integrator.run(n_timing_steps)
time_per_step = (time.time() - start) / float(n_timing_steps)
ns_per_hour = 3600.0 * time_step_fs * 1e-6 / time_per_step
ns_per_day = 24.0 * ns_per_hour
print("Yield:", ns_per_day, "ns/day")

if args.visu:
    visualizer.run(10)
else:
    print("\n-->Equilibration")
    n_int_steps = 10
    n_int_cycles = 100

    for i in range(n_int_cycles):
        system.integrator.run(n_int_steps)
        print("{0:.0f} %".format(100.0 * i / n_int_cycles))

    print("\n-->Integration")

    n_int_steps = 1000
    n_int_cycles = int(args.walltime * 3600.0 / time_per_step / n_int_steps)
    print(
        "Simulating for", args.walltime, "h, which is", n_int_cycles, "cycles x",
        n_int_steps, "steps, which is", args.walltime * ns_per_hour, "ns simulation time")

    file_traj = open(args.path + "traj.xyz", "w")
    n_parts_tot = len(system.part)

    for i in range(n_int_cycles):
        system.integrator.run(n_int_steps)

        # XYZ TRAJECTORY
        file_traj.write(str(n_parts_tot) + '\n')
        file_traj.write(
            "t(ns) = " + str(time_step_ns * i * n_int_steps) + '\n')
        for p in system.part:
            file_traj.write(
                shortTypes[p.type] + " " + ' '.join(map(str, p.pos_folded)) + '\n')

        # FOR RDFS
        system.analysis.append()

        print("{0:.1f} %".format(100.0 * i / n_int_cycles))

    # RDFs
    rdf_bins = 100
    r_min = 0.0
    r_max = system.box_l[0] / 2.0
    r, rdf_00 = system.analysis.rdf(rdf_type='<rdf>',
                                    type_list_a=[types["PF6"]],
                                    type_list_b=[types["PF6"]],
                                    r_min=r_min,
                                    r_max=r_max,
                                    r_bins=rdf_bins)

    r, rdf_11 = system.analysis.rdf(rdf_type='<rdf>',
                                    type_list_a=[types["BMIM_COM"]],
                                    type_list_b=[types["BMIM_COM"]],
                                    r_min=r_min,
                                    r_max=r_max,
                                    r_bins=rdf_bins)
    r, rdf_01 = system.analysis.rdf(rdf_type='<rdf>',
                                    type_list_a=[types["PF6"]],
                                    type_list_b=[types["BMIM_COM"]],
                                    r_min=r_min,
                                    r_max=r_max,
                                    r_bins=rdf_bins)

    rdf_fp = open(args.path + "rdf.dat", 'w')
    for i in range(rdf_bins):
        rdf_fp.write("%1.5e %1.5e %1.5e %1.5e\n" %
                     (r[i], rdf_00[i], rdf_11[i], rdf_01[i]))
    rdf_fp.close()

    file_traj.close()

    print("-->Done")
