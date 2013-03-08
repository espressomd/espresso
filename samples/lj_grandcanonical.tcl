#############################################################
#
# Hybrid MD/MC:
# Lennard-Jones system at constant chemical potential
#
# For explanation of the MC particle/insertion moves, compare
# Frenkel&Smit Molecular Simulation, pp. 126 ff.
#
# This script also includes tail corrections for energy and
# pressure
#
#############################################################
#
# Copyright (C) 2010,2012,2013 The ESPResSo project
# Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 
#   Max-Planck-Institute for Polymer Research, Theory Group
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


#############################################################
# parameters

# simulation box size
set box_l  6.3

# the chemical potential in kT as difference to log(volume/L^3)kT,
# where L is the thermal wave length. That is:
# chem_pot = mu/kT + log(V/L^3)
#
# for the ideal gas: mu = log(N L^3 / V)kT => chem_pot = log(N)
# => excess chemical potential mu_ex = kT(chem_pot - log(N))
# or chem_pot = mu_ex/kT + log(N)

set chem_pot 1

# initial number of particles, simply assume ideal gas
set initial_density [expr exp($chem_pot)*pow($box_l, -3)]
# but don't put too many particles...
if {$initial_density > 0.8} { set initial_density 0.8 }

# particle-particle interactions
set lj_eps   1
set lj_sigma 1
set lj_cut   2.5
set lj_shift [calc_lj_shift $lj_sigma $lj_cut]

# integration parameters
set time_step 0.01

# reduced temperature
set temperature 2.0
# friction of Langevin thermostat
set gamma 1.0

# number of time steps per cycle during warmup
set warm_steps  1000
# number of cycles during warmup
set warm_cycles   10

# number of time steps per cycle during sampling
# note that after each cycle, a number of particle
# exchanges is tried
set int_steps     10
# number of cycles during sampling
set int_cycles 50000
# number of sampling cycles to skip as equilibration for averaging
set equ_cycles  1000
# number of exchange tries during a cycle
set int_nr_exchanges 2

# name of data file
set filename "lj-mu_${chem_pot}.data"

set tcl_precision 6

set random_seed [pid]

# end of parameters
#############################################################

#############################################################
# subroutines
#############################################################

# get the average temperature per particle
proc current_temperature {} {
    if {[setmd n_part] > 0} {
	return [expr [analyze energy kinetic]/(1.5*[setmd n_part])]
    } {
	return 0
    }
}

# get the potential energy of the system. N is the number of LJ particles
# for which to do the tail correction
proc potential_energy {N} {
    global V
    global lj_eps lj_sigma lj_cut lj_shift

    # insert tail corrections
    if {$N > 0} {
	# potential energy
	set eng [expr [analyze energy total] - [analyze energy kinetic]]

	set dens [expr double($N)/$V]

	# tail correction
	set tail [expr $N*8./3.*[PI]*$dens*$lj_eps*pow($lj_sigma, 3)*\
		      (1./3.*pow($lj_sigma/$lj_cut, 9) - pow($lj_sigma/$lj_cut, 3))]

	# shift correction
	set shift [expr -$N*8./3.*[PI]*$dens*$lj_eps*pow($lj_cut, 3)*\
		       $lj_shift]

	return [expr $eng + $tail + $shift]
    } {
	return 0
    }
}

# get the pressure with tail corrections
proc virial_pressure {} {
    global V
    global lj_eps lj_sigma lj_cut

    # virial pressure
    set press [expr [analyze pressure total] - [analyze pressure ideal]]

    # tail correction
    set dens [expr [setmd n_part]/$V]
    set tail [expr 16./3.*[PI]*pow($dens, 2)*$lj_eps*pow($lj_sigma, 3)*\
		  (2./3.*pow($lj_sigma/$lj_cut, 9) - pow($lj_sigma/$lj_cut, 3))]

    return [expr $press + $tail]
}

proc create_particle {} {
    global box_l
    set px [expr rand()*$box_l]
    set py [expr rand()*$box_l]
    set pz [expr rand()*$box_l]

    # add particle after last known particle
    set id [expr [setmd max_part] + 1]
    part $id pos $px $py $pz v 0 0 0

    return $id
}

proc create_particles {n} {
    global box_l
    for {set cnt 0} {$cnt < $n} {incr cnt} {
	create_particle
    }
}

proc delete_particle {id} {
    if {$id == [setmd max_part]} {
	# last particle, just delete
	part $id delete
    } {
	# otherwise, shift down last particle
	# this avoids that the particle identities get excessive
	# and with that some tables
	eval part $id pos [part [setmd max_part] pr pos]
	part [setmd max_part] delete
    }
}

proc acceptance_probability {denergy} {
    global beta

    # limit energy difference; output should not be larger than 1
    if { $denergy < 0 } { set denergy 0 }
    return [expr exp(-$beta*$denergy)]
}

proc try_particle_insert {} {
    global temperature chem_pot box_l

    set Eold [potential_energy [setmd n_part]]
    
    # place new particle
    set id [create_particle]

    set Nnew [setmd n_part]
    set Enew [potential_energy [setmd n_part]]

    # Boltzmann weight
    set bf [acceptance_probability [expr -$temperature*($chem_pot - log($Nnew)) + $Enew - $Eold]]

    if {rand() < $bf} {
	# accepted
	return 1
    } {
	# rejected, delete particle again
	delete_particle $id
	return 0
    }
}

proc try_particle_delete {} {
    global temperature chem_pot box_l

    if {[setmd n_part] < 1} { return 0 }

    set Eold [potential_energy [setmd n_part]]
    set Nold [setmd n_part]

    # TRICK: we do not really delete the particle, just give it another
    # type, so that it does not interact. This is easier than removing it
    # and potentially restoring it, if the move was rejected.

    # randomly pick a particle, and "remove"
    set id [expr int(rand()*([setmd max_part] + 1))]
    part $id type 1

    # CAVEAT! Since we make the particle disappear only through its type,
    # it still counts to n_part, but we must not do the tail correction for it!
    set Enew [potential_energy [expr [setmd n_part] - 1]]

    # fix the particles type again
    part $id type 0

    # Boltzmann weight
    set bf [acceptance_probability [expr $temperature*($chem_pot - log($Nold)) + $Enew - $Eold]]

    if {rand() < $bf} {
	# accept, really remove particle
	delete_particle $id
	return 1
    } {
	# reject, keep things as is
	return 0
    }
}

#############################################################
# main
#############################################################

# set up simulation box
setmd box_l $box_l $box_l $box_l
set V [expr pow($box_l, 3)]

# set up integrator
setmd time_step $time_step
setmd skin 0.4
thermostat langevin $temperature $gamma
set beta [expr 1./$temperature]

# particles
set N [expr int($V*$initial_density)]
puts "creating $N particles in a box of volume $V (density = $initial_density)"

create_particles $N

# interaction
inter 0 0 lennard-jones $lj_eps $lj_sigma $lj_cut $lj_shift 0

# warmup
#############################################################

set cap 5
for {set cycle 0} {$cycle < $warm_cycles} {incr cycle} {
    set cap [expr $cap + 5]
    inter forcecap $cap

    puts -nonewline "\r[setmd time] T=[current_temperature] E_pot=[potential_energy [setmd n_part]]"
    flush stdout

    integrate $warm_steps
}

inter forcecap 0

puts "\n warmup done"

# main sampling
#############################################################

# acceptance statistics
set insert_tries 0
set insert_acc   0
set delete_tries 0
set delete_acc   0

# averages
set sum_P 0
set sum_N 0

set results [open $filename "w"]
puts $results "\# time N E_tot Pressure"

for {set cycle 0} {$cycle < $int_cycles} {incr cycle} {
    # propagate
    integrate $int_steps

    # try particle insertion/deletion
    for {set try 0} {$try < $int_nr_exchanges} {incr try} {
	# insert/delete with equal probability
	if {[expr rand()] < .5} {
	    incr insert_tries
	    if {[try_particle_insert]} {incr insert_acc }
	} {
	    incr delete_tries
	    if {[try_particle_delete]} {incr delete_acc }
	}
    }

    # console information
    puts -nonewline [format "\r%8.1f N=%4d T=%8.3f E_pot=%8.3f" \
			 [setmd time] [setmd n_part] [current_temperature] \
			 [potential_energy [setmd n_part]]]
    if { $insert_tries > 0 && $delete_tries > 0} {
	puts -nonewline [format " rate_insert=%4.3f rate_delete=%4.3f" \
			     [expr double($insert_acc)/$insert_tries] \
			     [expr double($delete_acc)/$delete_tries]]
    }
    flush stdout

    # record observables
    if {[setmd n_part] > 0} {
	set etot [expr [potential_energy [setmd n_part]] + [setmd n_part]*1.5*$temperature]
	set P    [expr [virial_pressure]                 + [setmd n_part]/$V*$temperature]
    } {
	set etot 0
	set P    0
    }

    puts $results [format "%8.2f %5d %8.3f %8.3f" \
		       [setmd time] [setmd n_part] $etot $P]
    flush $results

    # update averages
    if {$cycle >= $equ_cycles} {
	set sum_N [expr $sum_N + [setmd n_part]]
	set sum_P [expr $sum_P + $P]
    }
}
close $results

puts "\ndone"
set sampling_cycles [expr $int_cycles - $equ_cycles]
puts "\naverage \#particles after 1000 cycles was: [expr double($sum_N)/$sampling_cycles]"
puts "\naverage    pressure after 1000 cycles was: [expr double($sum_P)/$sampling_cycles]"
