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
import numpy

print("""
=======================================================
=                    lj_liquid.py                     =
=======================================================

Program Information:""")
print(code_info.features())

dev = "cpu"

# System parameters
#############################################################

# 10 000  Particles
box_l = 10.7437
density = 0.7

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
#es._espressoHandle.Tcl_Eval('thermostat langevin 1.0 1.0')
system.thermostat.set_langevin(kT=1.0, gamma=1.0)

# warmup integration (with capped LJ potential)
warm_steps = 100
warm_n_times = 30
# do the warmup until the particles have at least the distance min__dist
min_dist = 0.9

# integration
int_steps = 1000
int_n_times = 5


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

volume = box_l * box_l * box_l
n_part = int(volume * density)

for i in range(n_part):
    system.part.add(id=i, pos=numpy.random.random(3) * system.box_l)

analyze.distto(system, 0)

print("Simulate {} particles in a cubic simulation box {} at density {}."
      .format(n_part, box_l, density).strip())
print("Interactions:\n")
act_min_dist = analyze.mindist(system)
print("Start with minimal distance {}".format(act_min_dist))

system.max_num_cells = 2744

#############################################################
#  Warmup Integration                                       #
#############################################################

# open Observable file
obs_file = open("pylj_liquid.obs", "w")
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
    act_min_dist = analyze.mindist(system)
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
set_file = open("pylj_liquid.set", "w")
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
    linear_momentum = analyze.analyze_linear_momentum(system=system)
    print(linear_momentum)
    # print(analyze.calc_rh(system,0,3,5))
    # print(analyze.calc_rg(system,0,3,5))
    # print(analyze.calc_re(system,0,3,5))

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
end_file = open("pylj_liquid.end", "w")
end_file.write("{ time %f } \n { box_l %f }\n" % (system.time, box_l))
end_file.write("{ particles {id pos type} }")
for i in range(n_part):
    end_file.write("%s\n" % system.part[i].pos)
    # id & type not working yet

obs_file.close()
set_file.close()
end_file.close()
# es._espressoHandle.die()

# terminate program
print("\nFinished.")
