#
# Copyright (C) 2013,2014,2015,2016 The ESPResSo project
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
from espressomd import thermostat
from espressomd import code_info
from espressomd import analyze
from espressomd import integrate
from espressomd import electrostatics
import numpy

print("""
=======================================================
=                    debye_hueckel.py                 =
=======================================================

Program Information:""")
print(code_info.features())

dev = "cpu"

# Constants
#############################################################
N_A = 6.022e23
pi = 3.14159265359

# System parameters
#############################################################

box_l = 10
# Molar salt concentration
mol_dens = 0.1
# Number density of ions
num_dens = mol_dens * N_A
# Convert to MD units with lj_sig = 7.14 Angstrom
num_dens = num_dens * 3.64e-25

volume = box_l * box_l * box_l
n_part = int(volume * num_dens)

# Interaction parameters (repulsive Lennard Jones)
#############################################################

lj_eps = 1.0
lj_sig = 1.0
lj_cut = 1.12246
lj_cap = 20

# Integration parameters
#############################################################
system = espressomd.System()
system.time_step = 0.01
system.skin = 0.4

thermostat.Thermostat().set_langevin(1.0, 1.0)

# warmup integration (with capped LJ potential)
warm_steps = 100
warm_n_times = 30
# do the warmup until the particles have at least the distance min__dist
min_dist = 0.9

# integration
int_steps = 1000
int_n_times = 10


#############################################################
#  Setup System                                             #
#############################################################

# Interaction setup
#############################################################

system.box_l = [box_l, box_l, box_l]

system.non_bonded_inter[0, 0].lennard_jones.set_params(
    epsilon=lj_eps, sigma=lj_sig,
    cutoff=lj_cut, shift="auto")
system.non_bonded_inter.set_force_cap(lj_cap)


print("LJ-parameters:")
print(system.non_bonded_inter[0, 0].lennard_jones.get_params())

# Particle setup
#############################################################

for i in range(n_part):
    system.part.add(id=i, pos=numpy.random.random(3) * system.box_l)

for i in range(n_part / 2):
    system.part[2 * i].q = -1.0
    system.part[2 * i].type = 1
    system.part[2 * i + 1].q = 1.0
    system.part[2 * i + 1].type = 2

# for i in range(n_part-1):
#    print("Particle {} has charge {} and is of type {}.".format(i,system.part[i].q,system.part[i].type))

# Activating the Debye-Hueckel interaction
# The Bjerrum length is set to one. Assuming the solvent is water, this
# means that lj_sig is 0.714 nm in SI units.
l_B = 1
# inverse Debye length for 1:1 electrolyte in water at room temperature (nm)
dh_kappa = numpy.sqrt(mol_dens) / 0.304
# convert to MD units
dh_kappa = dh_kappa / 0.714
dh = electrostatics.CDH(bjerrum_length=l_B, kappa=dh_kappa,
                        r_cut=5 / dh_kappa, eps_int=1, eps_ext=1, r0=0.5, r1=1, alpha=1)
system.actors.add(dh)
print(system.actors)

analyze.distto(system, 0)

print("Simulate monovalent salt in a cubic simulation box {} at molar concentration {}."
      .format(box_l, mol_dens).strip())
print("Interactions:\n")
act_min_dist = analyze.mindist(es)
print("Start with minimal distance {}".format(act_min_dist))

system.max_num_cells = 2744

#############################################################
#  Warmup Integration                                       #
#############################################################

# open Observable file
obs_file = open("pydebye_hueckel.obs", "w")
obs_file.write("# Time\tE_tot\tE_kin\tE_pot\n")
# set obs_file [open "$name$ident.obs" "w"]
# puts $obs_file "\# System: $name$ident"
# puts $obs_file "\# Time\tE_tot\tE_kin\t..."

print("""
Start warmup integration:
At maximum {} times {} steps
Stop if minimal distance is larger than {}
""".strip().format(warm_n_times, warm_steps, min_dist))

# set LJ cap
lj_cap = 20
system.non_bonded_inter.set_force_cap(lj_cap)
print(system.non_bonded_inter[0, 0].lennard_jones)

# Warmup Integration Loop
i = 0
while (i < warm_n_times and act_min_dist < min_dist):
    integrate.integrate(warm_steps)
    # Warmup criterion
    act_min_dist = analyze.mindist(es)
#  print("\rrun %d at time=%f (LJ cap=%f) min dist = %f\r" % (i,system.time,lj_cap,act_min_dist), end=' ')
    i += 1

#   write observables
#    puts $obs_file "{ time [setmd time] } [analyze energy]"

#   Increase LJ cap
    lj_cap = lj_cap + 10
    system.non_bonded_inter.set_force_cap(lj_cap)

# Just to see what else we may get from the c code
print("""
ro variables:
cell_grid     {0.cell_grid}
cell_size     {0.cell_size} 
local_box_l   {0.local_box_l} 
max_cut       {0.max_cut}
max_part      {0.max_part}
max_range     {0.max_range} 
max_skin      {0.max_skin}
n_nodes       {0.n_nodes}
n_part        {0.n_part}
n_part_types  {0.n_part_types}
periodicity   {0.periodicity}
transfer_rate {0.transfer_rate}
verlet_reuse  {0.verlet_reuse}
""".format(system))

# write parameter file

# polyBlockWrite "$name$ident.set" {box_l time_step skin} ""
set_file = open("pydebye_hueckel.set", "w")
set_file.write("box_l %s\ntime_step %s\nskin %s\n" %
               (box_l, system.time_step, system.skin))

#############################################################
#      Integration                                          #
#############################################################
print("\nStart integration: run %d times %d steps" % (int_n_times, int_steps))

# remove force capping
lj_cap = 0
system.non_bonded_inter.set_force_cap(lj_cap)
print(system.non_bonded_inter[0, 0].lennard_jones)

# print initial energies
#energies = es._espressoHandle.Tcl_Eval('analyze energy')
energies = analyze.energy(system=system)
print(energies)

j = 0
for i in range(0, int_n_times):
    print("run %d at time=%f " % (i, system.time))

#  es._espressoHandle.Tcl_Eval('integrate %d' % int_steps)
    integrate.integrate(int_steps)

#  energies = es._espressoHandle.Tcl_Eval('analyze energy')
    energies = analyze.energy(system=system)
    print(energies)
    obs_file.write('{ time %s } %s\n' % (system.time, energies))

#   write observables
#    set energies [analyze energy]
#    puts $obs_file "{ time [setmd time] } $energies"
#    puts -nonewline "temp = [expr [lindex $energies 1 1]/(([degrees_of_freedom]/2.0)*[setmd n_part])]\r"
#    flush stdout

#   write intermediate configuration
#    if { $i%10==0 } {
#	polyBlockWrite "$name$ident.[format %04d $j]" {time box_l} {id pos type}
#	incr j
#    }

# write end configuration
end_file = open("pydebye_hueckel.end", "w")
end_file.write("{ time %f } \n { box_l %f }\n" % (system.time, box_l))
end_file.write("{ particles {type q pos} }")
for i in range(n_part - 1):
    end_file.write("%s\t%s\t%s\n" %
                   (system.part[i].type, system.part[i].q, system.part[i].pos))
    # id & type not working yet

obs_file.close()
set_file.close()
end_file.close()
# es._espressoHandle.die()

# terminate program
print("\nFinished.")
